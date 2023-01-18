#!/usr/bin/env sh
xvfb-run julia --project=@. -e 'using Pkg; Pkg.build(); Pkg.test(coverage=true)'