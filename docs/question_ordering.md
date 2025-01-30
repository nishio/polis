# Polisの質問提示順序の仕組み

## 概要
Polisでは、ユーザーに質問（コメント）を提示する順序は、単純なランダム表示ではなく、情報量の最大化を目指した複数の要因によって決定されています。

## コアロジック

### 優先度の計算（Priority Metric）
質問の優先度は`priority-metric`関数によって計算されます：

ファイル: `math/src/polismath/math/conversation.clj`
```clojure
;; メタ質問の優先度（固定値）
(def meta-priority 7)

(defn priority-metric
  [is-meta A P S E]
  (matrix/pow
    (if is-meta
      meta-priority  ; メタ質問は常に高優先度
      (* (importance-metric A P S E)
         ;; 新しいコメントを上位に表示するためのスケーリング係数
         (+ 1 (* 8 (matrix/pow 2 (/ S -5))))))
    2))
```

### 重要度の計算（Importance Metric）
各質問の基本的な重要度は以下の要因を考慮して計算されます：

ファイル: `math/src/polismath/math/conversation.clj`
```clojure
(defn importance-metric
  [A P S E]
  (let [p (/ (+ P 1) (+ S 2))
        a (/ (+ A 1) (+ S 2))]
    (* (- 1 p) (+ E 1) a)))
```

パラメータの説明：
- `A`: 賛成票の数（Agrees）
- `P`: パス（スキップ）の数（Passes）
- `S`: 総投票数（Total Votes）
- `E`: 意見の極端さ（Extremity）

### 優先順位の決定要因

1. **メタ質問の優先度**
   - メタ質問（`is-meta = true`）は常に高い固定優先度（7）が与えられます

2. **新規性の考慮**
   - 新しい質問は指数関数的なスケーリング係数により優先的に表示されます
   - スケーリング係数: `(+ 1 (* 8 (matrix/pow 2 (/ S -5))))`
   - 投票数（S）が少ないほど、この係数は大きくなります

3. **投票パターンの影響**
   - 賛成票（A）が多いコメントは優先度が上がります
   - パス/スキップ（P）が多いコメントは優先度が下がります
   - これらは`importance-metric`関数内で計算されます

4. **意見の極端さ**
   - より極端な意見（E: extremity）を持つコメントは優先される傾向にあります
   - extremityは意見空間でのコメントの位置に基づいて計算されます

## 質問提示順序の詳細な仕組み

### 1. Extremityの計算方法
ファイル: `math/src/polismath/math/pca.clj`
```clojure
(defn with-proj-and-extremtiy [pca]
  (let [cmnt-proj (pca-project-cmnts pca)
        cmnt-extremity
        (mapv (fn [row] (matrix/length row))
          (matrix/rows cmnt-proj))]
    (assoc pca
           :comment-projection (matrix/transpose cmnt-proj)
           :comment-extremity cmnt-extremity)))
```

Extremity（意見の極端さ）は以下のように計算されます：

1. **PCAによる次元削減**
   - 各コメントを低次元空間（主成分空間）に射影
   - この空間では「座標の原点付近」が意見分布の中心・平均的な位置を表す

2. **距離による極端さの定量化**
   - 射影された各コメントの座標（ベクトル）のユークリッド距離（L2ノルム）を計算
   - この距離が「平均的立ち位置からどれだけ離れているか＝極端さ」を表す

3. **革新的な特徴**
   - 複数の主成分全体での離れ具合を一本化して評価可能
   - 「どの方向に極端か」よりも「どれだけ極端か」という量的評価に重点
   - 高次元の意見データでも比較的シンプルに扱える

### 2. メタ質問の定義と優先度
メタ質問（`is-meta = true`）は以下のような特徴を持ちます：

1. **定義**
   - システム自体の運用や会話の構造に関する"メタ的な話題"
   - 運用側の設定やユーザーフラグ（`:meta-tids`への登録）で判定

2. **固定優先度7の理由**
   - システム全体の設計・運営に関わる重要な話題として特別扱い
   - 通常の優先度計算から独立させることで、一定の注目度を確保
   - サポートやモデレーション方針など、利用者全体の利便性に関わる議論を保護

3. **設計上の工夫**
   - 通常の質問とは別枠で扱うことで、埋もれや過度な優先を防止
   - 運用ポリシーとしての"人為的優先度"をシステムに組み込む

### 3. 新規性係数のチューニング
```clojure
(+ 1 (* 8 (matrix/pow 2 (/ S -5))))
```

この式における定数は以下の意図で設定されています：

1. **投票数による減衰（-5）**
   - 5票増えるごとにスケーリング係数が約1/2になる半減期的な挙動
   - 投票数（S）の増加に応じて指数関数的に影響を減少

2. **新規コメントのブースト係数（8）**
   - 新規コメントを一時的に大きくブーストする係数
   - 実運用テストを通じて調整された経験的な値

3. **設計の特徴**
   - 指数関数による滑らかな制御
   - 新規投稿の可視性を確保しつつ、継続的に支持されるコメントとのバランスを取る
   - 運用者が感覚的に調整しやすい数値設定

## 関連コード

### コメント優先度の計算（Comment Priorities）
ファイル: `math/src/polismath/math/conversation.clj`
```clojure
:comment-priorities
(plmb/fnk [conv group-votes pca tids meta-tids]
  (let [group-votes (:group-votes conv)
        extremities (into {} (map vector tids (:comment-extremity pca)))]
    (plmb/map-from-keys
      (fn [tid]
        (let [{:as total-votes :keys [A D S P]}
              (reduce
                (fn [votes [gid data]]
                  (let [{:as data :keys [A S D] :or {A 0 S 0 D 0}} (get-in data [:votes tid])
                        data (assoc data :P (+ (- S (+ A D))))]
                    (reduce
                      (fn [votes' [k v]]
                        (update votes' k + v))
                      votes
                      data)))
                {:A 0 :D 0 :S 0 :P 0}
                group-votes)
              extremity (or (get extremities tid)
                            (do
                              (log/warn "No extremity for tid" tid "zid" (:zid conv))
                              0))]
          (priority-metric (meta-tids tid) A P S extremity)))
      tids)))
```

このコードは、各コメントの優先度を計算し、それに基づいて質問の提示順序を決定します。
