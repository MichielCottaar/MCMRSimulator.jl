#!/usr/bin/env sh
# Install xvfb
apt-get update
apt-get install -y xvfb

# run tests
xvfb-run julia --project=@. -e 'using Pkg; Pkg.build(); Pkg.test(coverage=true)'