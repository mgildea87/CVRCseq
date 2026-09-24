FROM continuumio/miniconda3:latest

COPY workflow/envs/CVRCseq.yml /env.yml

RUN set -eux; \
		conda env create -f /env.yml -n CVRCseq; \
		conda clean -a -y

ENV PATH=/opt/conda/envs/CVRCseq/bin:$PATH
ENV CONDA_DEFAULT_ENV=CVRCseq

# --- PATCH THE SHELL SCRIPT HERE ---
# Copy your fixed .sh file over the package's version
COPY patches/SEACR_1.3.sh /opt/conda/envs/CVRCseq/bin/SEACR_1.3.sh

# Ensure the script has proper execution permissions
RUN chmod +x /opt/conda/envs/CVRCseq/bin/SEACR_1.3.sh