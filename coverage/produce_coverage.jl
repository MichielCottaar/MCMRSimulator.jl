#!/usr/bin/env julia --project=@coverage
using Coverage
coverage = process_folder()
c, t = get_summary(coverage)
using Printf
@printf "Test coverage %.2f%%\n" 100c / t

open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end