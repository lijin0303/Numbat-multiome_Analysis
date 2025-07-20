# Updated epiAneufinder Docker image
# Based on existing scatac_tool but with updated epiAneufinder package

FROM gcr.io/broad-getzlab-workflows/scatac_tool:v41

LABEL maintainer="Ruitong Li <ruitong@broadinstitute.org>"
LABEL description="Updated epiAneufinder package from GitHub commit 6969d6c"

# Update system packages
RUN apt-get update && apt-get install -y \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install/update R packages needed for the latest epiAneufinder
RUN R -e "install.packages(c('devtools', 'remotes'), repos='https://cloud.r-project.org')"

# Remove old epiAneufinder and install the specific GitHub version
RUN R -e "remove.packages('epiAneufinder', lib=.libPaths()[1])" || true

# Install the specific commit of epiAneufinder
RUN R -e "devtools::install_github('colomemaria/epiAneufinder', ref='6969d6c4926acb707087f7a03dcfabca7edd5cfe', force=TRUE, dependencies=TRUE)"

# Verify installation
RUN R -e "library(epiAneufinder); packageVersion('epiAneufinder')"

# Set working directory
WORKDIR /workspace

# Default command
CMD ["/bin/bash"]
