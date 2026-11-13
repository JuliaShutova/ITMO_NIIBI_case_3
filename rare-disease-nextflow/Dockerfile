FROM mambaorg/micromamba:1.5.10
COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n pipeline -f /tmp/environment.yml && micromamba clean --all --yes
ENV PATH=/opt/conda/envs/pipeline/bin:$PATH
WORKDIR /work
