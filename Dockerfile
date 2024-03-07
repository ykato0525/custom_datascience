FROM jupyter/datascience-notebook
WORKDIR /work
EXPOSE 8888
COPY . ./

# jupyter lab extensions.
# RUN conda update conda -y
RUN conda install -c conda-forge scanpy python-igraph leidenalg
RUN conda install -c bioconda gseapy
# RUN conda install -c conda-forge numpy==1.23.1
RUN mamba install -c bioconda star=2.7.5a -y
RUN mamba install -c bioconda fastqc -y
RUN mamba install -c bioconda fastp -y
RUN mamba install -c bioconda samtools -y
RUN mamba install -c bioconda rsem -y
RUN mamba install -c bioconda multiqc -y
RUN mamba install -c bioconda bowtie -y
#scanpyとのバージョンを合わせるため
RUN pip install nummpy==1.23.1
# NGS解析をjupyter上で可視化するために入れておきました。
RUN pip install igv-notebook -y
# 以下ChIP-seqのためのライブラリです
RUN conda conda create -n chip python=3.9.6
RUN conda activate chip
RUN mamba install -c bioconda bowtie -y
RUN mamba install -c bioconda MACS2 -y
RUN conda deactivate
RUN Rscript /work/setup.R
RUN bash setup.sh
