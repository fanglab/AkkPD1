
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
# Set locale to prevent LC_* warnings
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && apt-get install -y \
    software-properties-common dirmngr gnupg curl ca-certificates \
    libcurl4-openssl-dev libssl-dev libxml2-dev libxt-dev \
    libjpeg-dev libpng-dev zlib1g-dev build-essential \
    && curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /etc/apt/trusted.gpg.d/cran.gpg \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    && apt-get update && apt-get install -y r-base

RUN R -e "install.packages(c('caret', 'readr', 'dplyr', 'MLeval', 'optparse', 'pROC', 'doParallel', 'randomForest'), repos='http://cran.rstudio.com/')"

RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='http://cran.rstudio.com/'); BiocManager::install('MLeval')"

COPY akk_gene_ML/rf_model_training_pipeline.R /app/rf_model_training_pipeline.R
WORKDIR /app
ENTRYPOINT ["Rscript", "rf_model_training_pipeline.R"]
