FROM julia:latest
WORKDIR /env

COPY . .

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=. -e 'import MCMRSimulator.CLI: run_main; run_main(["geometry", "create", "walls"])'; exit 0
RUN julia --project=. -e 'import MCMRSimulator.CLI: run_main; run_main(["run"])'; exit 0
RUN julia --project=. -e 'import MCMRSimulator.CLI: run_main; run_main(["sequence", "dwi"])'; exit 0
ENTRYPOINT [ "julia", "--project=/env", "-e", "import MCMRSimulator.CLI: run_main; run_main()", "--" ]