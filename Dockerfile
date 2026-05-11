FROM continuumio/miniconda3:latest

COPY workflow/envs/CVRCseq.yml /env.yml

RUN set -eux; \
		conda env create -f /env.yml -n CVRCseq; \
		seacr_sh="$(find /opt/conda/envs/CVRCseq -type f -name 'SEACR_1.3.sh' | head -n 1)"; \
		if [ -z "$seacr_sh" ]; then \
			echo 'SEACR_1.3.sh not found in CVRCseq environment' >&2; \
			exit 1; \
		fi; \
		perl -0pi -e 's#^path=`dirname \$0`#path="\$(cd "\$(dirname "\$(command -v "\$0")")" && pwd)"#m; s#\./SEACR_1\.3\.R#\$path/SEACR_1.3.R#g; s#(?<![/[:alnum:]_.-])SEACR_1\.3\.R#\$path/SEACR_1.3.R#g' "$seacr_sh"; \
		conda clean -a -y

ENV PATH=/opt/conda/envs/CVRCseq/bin:$PATH
ENV CONDA_DEFAULT_ENV=CVRCseq
