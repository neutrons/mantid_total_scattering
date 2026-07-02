# Start from a clean conda-forge base image with mamba installed
FROM condaforge/miniforge3:latest

# Install system dependencies required for Mantid and other tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    git \
    libgl1 \
    libglu1-mesa \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install pixi
RUN curl -fsSL https://pixi.sh/install.sh | sh

# Add pixi to PATH
ENV PATH="/root/.pixi/bin:${PATH}"

# Set working directory
WORKDIR /app

# Copy the project first (before pixi install)
COPY . .

# Install dependencies using pixi
# This will use the pixi.lock file if present, and install all features
# including mantidworkbench (>=6.14) and test dependencies
RUN pixi install --locked

# Activate pixi environment for all subsequent operations
# The pixi environment is already in PATH after 'pixi install'
ENV PATH="/app/.pixi/envs/default/bin:${PATH}"

# Set environment variables for mantid total scattering
ENV PYTHONPATH="/app:${PYTHONPATH}"

# Define the default command to run pytest
CMD ["pixi", "run", "mantidtotalscattering", "--help"]
