#!/usr/bin/env sh
apt-get update
apt-get install -y git
|
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'
mkdir -p public
mv docs/build public/$1