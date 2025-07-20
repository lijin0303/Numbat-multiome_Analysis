#!/bin/bash

# Build and push updated epiAneufinder Docker image

set -euo pipefail

# Configuration
IMAGE_NAME="gcr.io/broad-getzlab-workflows/scatac_tool"
NEW_TAG="v42-epianeu-updated"
DOCKERFILE="Dockerfile"

echo "=== Building updated epiAneufinder Docker image ==="
echo "Image: ${IMAGE_NAME}:${NEW_TAG}"
echo "Dockerfile: ${DOCKERFILE}"
echo

# Build the image
echo "Building Docker image..."
docker build -t "${IMAGE_NAME}:${NEW_TAG}" -f "${DOCKERFILE}" .

if [ $? -eq 0 ]; then
    echo "✓ Docker image built successfully"
else
    echo "✗ Docker build failed"
    exit 1
fi

# Test the image
echo
echo "Testing the new image..."
docker run --rm "${IMAGE_NAME}:${NEW_TAG}" R -e "
library(epiAneufinder)
cat('epiAneufinder version:', as.character(packageVersion('epiAneufinder')), '\n')
cat('Testing epiAneufinder function availability...\n')
if (exists('epiAneufinder', where='package:epiAneufinder')) {
    cat('✓ epiAneufinder function is available\n')
} else {
    cat('✗ epiAneufinder function not found\n')
    quit(status=1)
}
"

if [ $? -eq 0 ]; then
    echo "✓ Image test passed"
else
    echo "✗ Image test failed"
    exit 1
fi

# Push the image
echo
echo "Pushing image to registry..."
docker push "${IMAGE_NAME}:${NEW_TAG}"

if [ $? -eq 0 ]; then
    echo "✓ Image pushed successfully"
    echo
    echo "New image available at: ${IMAGE_NAME}:${NEW_TAG}"
    echo
    echo "To use the new image, update your scripts to use:"
    echo "  ${IMAGE_NAME}:${NEW_TAG}"
else
    echo "✗ Image push failed"
    exit 1
fi

echo
echo "=== Build and push completed successfully ==="
