#!/usr/bin/env zsh
julia --project=@. -e 'using Pkg; Pkg.test(coverage=true)'
julia --project=coverage coverage/produce_coverage.jl
rm **/*.cov
