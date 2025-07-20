#!/bin/bash
# Usage: bash run_congasp_docker.sh
# This script launches the Docker container, mounts $PWD/congasp as /workspace, and runs run_congasp_training.sh inside the container.
# 
# Before running, replace 'yourusername' with your actual Docker Hub username:
# DOCKER_IMAGE="yourusername/congas:latest"


DOCKER_IMAGE="ruitongl/congasp:latest"

# Run the script inside the container with increased memory limits
docker run --rm \
  --memory=16g \
  --memory-swap=32g \
  -v "$PWD/congasp:/workspace" \
  $DOCKER_IMAGE \
  bash /workspace/run_congasp_training.sh 