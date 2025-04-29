# Start from Ubuntu
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH

# Install basic Linux packages
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    bzip2 \
    git \
    build-essential \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda manually
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -t -i -p -y

# Set working directory
WORKDIR /pipeline
te environment for all future RUN/CMD
SHELL ["conda", "run", "-n", "akk_pipeline", "/bin/bash", "-c"]
# Copy environment.yml
COPY environment.yml .

# Create conda environment
RUN conda env create -f environment.yml

# Activate environment for all future RUN/CMD
SHELL ["conda", "run", "-n", "akk_pipeline", "/bin/bash", "-c"]

# Create data and DB directories
RUN mkdir -p /pipeline/data /pipeline/DB

# Copy pipeline scripts into container
COPY run_akk_pipeline.sh .
COPY assign_phylogroups.R .
COPY run_pcoa_plot.R .
COPY run_tests.sh .

# Make scripts executable
RUN chmod +x run_akk_pipeline.sh run_tests.sh

# Reset shell for runtime
SHELL ["/bin/bash", "-c"]

# Activate conda env and open bash when container runs
ENTRYPOINT ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate akk_pipeline && exec bash"]
