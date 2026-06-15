FROM mambaorg/micromamba:2.8.1 AS builder

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

RUN micromamba install -y -n base \
      -c conda-forge \
      -c bioconda \
      --override-channels \
      --strict-channel-priority \
      python=3.13 \
      pip \
      hatchling \
      git \
      blast \
      kma \
      "cgecore>=2,<3" \
      "biopython>=1.81" \
      "tabulate>=0.8" \
      platformdirs \
    && micromamba clean -a -y \
    && rm -rf /opt/conda/pkgs /root/.mamba/pkgs

ENV PATH=/opt/conda/bin:$PATH

WORKDIR /tmp/pmlst-src

COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN python -m pip install . --no-deps --no-build-isolation \
    && pmlst-download-db /opt/pmlst/db \
    && rm -rf /opt/pmlst/db/.git /tmp/pmlst-src /root/.cache/pip

FROM mambaorg/micromamba:2.8.1

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

RUN micromamba install -y -n base \
      -c conda-forge \
      -c bioconda \
      --override-channels \
      --strict-channel-priority \
      python=3.13 \
      blast \
      kma \
      "cgecore>=2,<3" \
      "biopython>=1.81" \
      "tabulate>=0.8" \
      platformdirs \
    && micromamba clean -a -y \
    && rm -rf /opt/conda/pkgs /root/.mamba/pkgs

ENV PATH=/opt/conda/bin:$PATH
ENV PMLST_DB=/opt/pmlst/db

COPY --from=builder /opt/conda/bin/pmlst /opt/conda/bin/pmlst
COPY --from=builder /opt/conda/bin/pmlst.py /opt/conda/bin/pmlst.py
COPY --from=builder /opt/conda/bin/pmlst-download-db /opt/conda/bin/pmlst-download-db
COPY --from=builder /opt/conda/lib/python3.13/site-packages/cli.py /opt/conda/lib/python3.13/site-packages/cli.py
COPY --from=builder /opt/conda/lib/python3.13/site-packages/pmlst /opt/conda/lib/python3.13/site-packages/pmlst
COPY --from=builder /opt/conda/lib/python3.13/site-packages/pmlst-2.1.0.dist-info /opt/conda/lib/python3.13/site-packages/pmlst-2.1.0.dist-info
COPY --from=builder /opt/pmlst/db /opt/pmlst/db

WORKDIR /data

CMD ["pmlst", "--help"]
