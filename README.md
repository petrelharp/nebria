- run SLiM in the `small/` directory to see the simulation on a smaller chunk of the landscape

- run `Rscript make_maps.R` to make the maps that SLiM uses

- `.tif` files: Yi-Ming says (3 June 2021) " The file names contain the name of the time period and the chunk of time. For example, BA13_15 means BA period that is from 15 kyr to 13 kyr."

- To run a simulation copy the required map files into a subdirectory,
    make a `params.json` JSON file with parameters that change the default,
    and run `slim ../nebria.slim` from within that directory.
