# hash:sha256:4fd3a0c690ec8f3715e38e25b36150ef7d008a4fe2689a00884ecaa091b4dcbd
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ghostscript=9.26~dfsg+0-0ubuntu0.18.04.12 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN /postInstall
