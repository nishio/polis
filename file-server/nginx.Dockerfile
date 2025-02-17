FROM docker.io/nginx:1.21.5-alpine

COPY nginx/nginx-ssl.site.default.conf /etc/nginx/conf.d/default.conf

# We only use these in testing.
COPY nginx/certs/snakeoil.cert.pem /etc/nginx/certs/snakeoil.cert.pem
COPY nginx/certs/snakeoil.key.pem  /etc/nginx/certs/snakeoil.key.pem

COPY nginx/certs/nhiro-fullchain.pem  /etc/nginx/certs/nhiro-fullchain.pem
COPY nginx/certs/nhiro-privkey.pem /etc/nginx/certs/nhiro-privkey.pem

EXPOSE 80
EXPOSE 443
