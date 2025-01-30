# Polisのクラスタリング実装の解説

## 概要
Polisでは、ユーザーの意見をグループ化するために2段階のk-meansクラスタリングを採用しています：

1. 最初に100個のベースクラスタを作成（base-clusters）
2. 次にそれらを2-5個のグループにまとめる（group-clusters）
3. シルエット係数を用いて最適なクラスタ数を選択

この文書では、以下の3つの観点からクラスタリングの実装を解説します：

1. クラスタリングの具体的な実装方法
2. データ追加時の計算コスト最適化手法
3. 階層的クラスタリングへの置き換えの可能性

## コアとなる実装

### 1. ベースクラスタリング（k=100）

ファイル: `/home/ubuntu/repos/polis/math/src/polismath/math/conversation.clj:427-450`
```clojure
:base-clusters
(plmb/fnk [conv proj-nmat in-conv opts']
  (let [in-conv-mat (nm/rowname-subset proj-nmat in-conv)]
    (sort-by :id
      (clusters/kmeans in-conv-mat
        (:base-k opts')           ; デフォルトで100
        :last-clusters (:base-clusters conv)
        :max-iters (:base-iters opts')))))
```

### 2. グループクラスタリング（k=2~5）

ファイル: `/home/ubuntu/repos/polis/math/src/polismath/math/conversation.clj:430-470`
```clojure
:group-clusterings
(plmb/fnk [conv base-clusters-weights base-clusters-proj opts']
    (plmb/map-from-keys
      (fn [k]
        (sort-by :id
          (clusters/kmeans base-clusters-proj k
            :last-clusters
              (when-let [last-clusterings (:group-clusterings conv)]
                (last-clusterings k))
            :cluster-iters (:group-iters opts')
            :weights base-clusters-weights)))
      (range 2 (inc (max-k-fn base-clusters-proj (:max-k opts'))))))
```

### 3. シルエット係数による評価

シルエット係数は、クラスタリングの品質を評価するための指標です。各データポイントについて、以下の計算を行います：

1. a(i): データポイントiと同じクラスタ内の他のポイントとの平均距離
2. b(i): データポイントiと最も近い他のクラスタの全ポイントとの平均距離
3. s(i) = (b(i) - a(i)) / max(a(i), b(i))

シルエット係数が高いほど、クラスタリングの品質が良いことを示します。

ファイル: `/home/ubuntu/repos/polis/math/src/polismath/math/clusters.clj:380-400`
```clojure
(defn silhouette
  "Compute the silhoette coefficient for either a cluster member, or for an entire clustering."
  ([distmat clusters member]
   (let [dist-row (nm/rowname-subset distmat [member])
         [a b]
         (reduce
           (fn [[a b] clst]
             (let [memb-clst? (some #{member} (:members clst))
                   membs (remove #{member} (:members clst))]
               (if (and memb-clst? (empty? membs))
                 (reduced [1 1])
                 (as-> membs data
                   (nm/colname-subset dist-row data)
                   (nm/get-matrix data)
                   (first data)
                   (matrix-stats/mean data)
                   (if memb-clst?
                     [data b]
                     [a (min data (or b data))])))))
           [nil nil]
           clusters)]
     (/ (- b a) (max b a))))
  ([distmat clusters]
   (let [cluster-distmat (nm/rowname-subset distmat (mapcat :members clusters))]
     (weighted-mean
       (map
         (partial silhouette distmat clusters)
         (nm/rownames cluster-distmat))))))

ファイル: `/home/ubuntu/repos/polis/math/src/polismath/math/conversation.clj:470-490`
```clojure
:group-clusterings-silhouettes
(plmb/fnk [group-clusterings bucket-dists]
  (plmb/map-vals (partial clusters/silhouette bucket-dists) group-clusterings))

:group-k-smoother
(plmb/fnk
  [conv group-clusterings group-clusterings-silhouettes opts']
  (let [{:keys [last-k last-k-count smoothed-k] :or {last-k-count 0}}
        (:group-k-smoother conv)
        count-buffer (:group-k-buffer opts')
        this-k       (apply max-key group-clusterings-silhouettes (keys group-clusterings))
        same         (if last-k (= this-k last-k) false)
        this-k-count (if same (+ last-k-count 1) 1)
        smoothed-k   (if (>= this-k-count count-buffer)
                       this-k
                       (if smoothed-k smoothed-k this-k))]
    {:last-k       this-k
     :last-k-count this-k-count
     :smoothed-k   smoothed-k}))
```

### 4. データ更新時の最適化

データが追加された際の計算コストを削減するため、以下の3つの主要な最適化手法が実装されています：

1. **クラスタの初期値最適化**
   - `clean-start-clusters`関数を使用して既存のクラスタ構造を維持
   - 新しいデータに対して既存のクラスタを再中心化
   - 必要な場合のみ新しいクラスタを追加

2. **部分的PCA更新**
   - ミニバッチ方式でPCAを更新
   - 学習率とフォーゲット率を用いて古い結果と新しい結果を統合
   - メモリ使用量と計算時間を削減

3. **サンプリングベースの更新**
   - データサイズに応じて100から1500の間でサンプルサイズを調整
   - ミニバッチPCAとの組み合わせで効率的な更新を実現

ファイル: `/home/ubuntu/repos/polis/math/src/polismath/math/clusters.clj:300-350`
```clojure
(defn clean-start-clusters
  "This function takes care of some possible messy situations which can crop up with using 'last-clusters'
  in kmeans computation, and generally gets the last set of clusters ready as the basis for a new round of
  clustering given the latest set of data."
  [data clusters k & {:keys [weights]}]
  (if (seq clusters)
    (let [clusters (safe-recenter-clusters data clusters :weights weights)
          uniq-clusters (uniqify-clusters clusters)
          possible-clusters (min k (count (distinct (into [] (matrix/rows (nm/get-matrix data))))))]
      (loop [clusters uniq-clusters]
        (let [clusters (recenter-clusters data clusters :weights weights)]
          (if (> possible-clusters (count clusters))
            (let [outlier (most-distal data clusters)]
              (if (> (:dist outlier) 0)
                (recur
                  (->
                    (mapv
                      (fn [clst]
                        (assoc clst :members
                          (remove (set [(:id outlier)]) (:members clst))))
                      clusters)
                    (conj {:id (inc (apply max (map :id clusters)))
                           :members [(:id outlier)]
                           :center (nm/get-row-by-name data (:id outlier))})))
                clusters))
            clusters))))
    (do
      (log/warn "Had to initialize clusters from clean-start-clusters.")
      (init-clusters data k))))

ファイル: `/home/ubuntu/repos/polis/math/src/polismath/math/conversation.clj:600-630`
```clojure
(defn partial-pca
  [mat pca indices & {:keys [n-comps iters learning-rate]
                      :or {n-comps 2 iters 10 learning-rate 0.01}}]
  (let [rating-subset (utils/filter-by-index mat indices)
        part-pca (pca/powerit-pca rating-subset n-comps
                     :start-vectors (:comps pca)
                     :iters iters)
        forget-rate (- 1 learning-rate)
        learn (fn [old-val new-val]
                (let [old-val (matrix/join old-val (repeat (- (matrix/dimension-count new-val 0)
                                                              (matrix/dimension-count old-val 0)) 0))]
                  (+ (* forget-rate old-val) (* learning-rate new-val))))]
    (fn [pca']
      {:center (learn (:center pca') (:center part-pca))
       :comps  (mapv #(learn %1 %2) (:comps pca') (:comps part-pca))})))
```

## Q1: クラスタリングの実装詳細

クラスタリングは以下の3段階で実行されます：

1. **ベースクラスタリング（k=100）**
   - `base-clusters`関数で実行
   - PCAで次元削減したデータを入力として使用
   - 前回のクラスタリング結果（`:last-clusters`）を初期値として利用可能
   - デフォルトで最大100回のイテレーション

2. **グループクラスタリング（k=2~5）**
   - `group-clusterings`関数で実行
   - ベースクラスタの重み付き中心を入力として使用
   - k=2から最大5（または`max-k-fn`による計算値）までの各kについてクラスタリングを実行
   - 各kについて前回の結果を初期値として利用可能

3. **最適なkの選択**
   - シルエット係数を計算（`group-clusterings-silhouettes`）
   - 最も高いシルエット係数を持つkを選択
   - 急激な変化を防ぐため`group-k-smoother`でスムージング
   - バッファ回数（デフォルト4回）連続で同じkが選ばれるまで前回のkを維持

## Q2: データ追加時の計算コスト最適化

システムは以下の3つの主要な最適化手法を採用しています：

1. **クラスタの初期値最適化**
   - 前回のクラスタリング結果を初期値として使用
   - `clean-start-clusters`関数で古いクラスタを新しいデータに適応
   - これにより収束までのイテレーション回数を削減

2. **部分的PCA更新**
   - `partial-pca`関数でミニバッチ更新を実装
   - 新しいデータの一部のみを使用してPCAを更新
   - 学習率とフォーゲット率を用いて古い結果と新しい結果を統合
   - メモリ使用量と計算時間を削減

3. **サンプリングベースの更新**
   - 大規模な会話では`sample-size-fn`でサンプルサイズを決定
   - データサイズに応じて100から1500の間でサンプルサイズを調整
   - ミニバッチPCAとの組み合わせで効率的な更新を実現

## Q3: 階層的クラスタリングへの置き換え可能性

現在のk-meansを階層的クラスタリングに置き換えることは技術的には可能ですが、以下の点を考慮する必要があります：

1. **インターフェースの互換性**
   - 現在のシステムは`:id`、`:members`、`:center`を持つクラスタ構造を前提
   - 階層的クラスタリングでもこの構造を維持する必要あり
   - `xy-clusters-to-nmat`などの関数は修正が必要

2. **シルエット係数の計算**
   - 現在のシルエット係数計算は任意のクラスタリングに対応
   - 階層的クラスタリングでも同じ評価方法が利用可能
   - ただし、デンドログラムのカット方法の検討が必要

3. **増分更新への対応**
   - 現在の`clean-start-clusters`のような最適化が使えない
   - 階層的クラスタリング用の増分更新アルゴリズムの実装が必要
   - 計算コストが増加する可能性が高い

4. **実装上の課題**
   - 現在の2段階クラスタリングを単一の階層的クラスタリングに置き換え
   - クラスタ数の自動決定方法の再設計が必要
   - `group-k-smoother`に相当する安定化機構の実装が必要

結論として、置き換えは可能ですが、以下の理由から慎重な検討が必要です：
- 既存コードベースへの大規模な変更が必要
- 増分更新の効率性が低下する可能性
- クラスタ数の決定方法の再設計が必要
- 現在のk-meansベースの実装は既に十分に最適化されている
