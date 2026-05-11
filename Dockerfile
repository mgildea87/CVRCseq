FROM continuumio/miniconda3:latest

COPY workflow/envs/CVRCseq.yml /env.yml

RUN set -eux; \
		conda env create -f /env.yml -n CVRCseq; \
		conda clean -a -y

ENV PATH=/opt/conda/envs/CVRCseq/bin:$PATH
ENV CONDA_DEFAULT_ENV=CVRCseq