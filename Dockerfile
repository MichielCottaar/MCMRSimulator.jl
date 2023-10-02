FROM julia:latest
WORKDIR /env

COPY . .

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=. -e 'import MCMRSimulator.CLI: run_main; run_main("geometry"); run_main("run"); run_main("sequence")'
ENTRYPOINT [ "julia", "--project=/env", "-e", "import MCMRSimulator.CLI: run_main; run_main()", "--" ]