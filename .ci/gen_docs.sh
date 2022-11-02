#!/usr/bin/env sh
apt-get update
apt-get install -y git

julia --project=docs docs/make.jl
mkdir -p public
mv docs/build public/$1