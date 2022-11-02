#!/usr/bin/env sh
julia -e 'using Pkg; Pkg.add("Coverage"); using Coverage; c, t = get_summary(process_folder()); using Printf; @printf "Test coverage %.2f%%\n" 100c / t'
