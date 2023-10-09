#!/usr/bin/env sh
set -e
tag=$1

function docker_tag_exists() {
    curl --silent -f -lSL https://hub.docker.com/v2/repositories/michielc/mcmr/tags/$1 > /dev/null 
}

if docker_tag_exists $tag ; then
    echo "Image with tag '$tag' already exists on docker hub."
    exit 1
fi
if docker_tag_exists $tag-amd64 ; then
    echo "Checking out tag '$tag' from git"
    git checkout $tag
    echo "Building $tag-arm64 locally"
    docker login
    docker build -t michielc/mcmr:$tag-arm64 .
    docker push michielc/mcmr:$tag-arm64
    docker manifest create michielc/mcmr:$tag michielc/mcmr:$tag-arm64 michielc/mcmr:$tag-arm64
    docker manifest push michielc/mcmr:$tag
else
    echo "ARM64 image for '$tag' has not been uploaded yet. Please wait for the gitlab pipeline to run before running this script."
    exit 1
fi


