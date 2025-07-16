#!/bin/bash
# Usage: bash run_congasp_docker.sh
# This script launches the Docker container, mounts $PWD/congasp as /workspace, and runs run_congasp_training.sh inside the container.

DOCKER_IMAGE="congasp:latest"

# Run the script inside the container
docker run --rm \
  -v "$PWD/congasp:/workspace" \
  $DOCKER_IMAGE \
  bash /workspace/run_congasp_training.sh 