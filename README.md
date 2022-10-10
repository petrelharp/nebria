- run SLiM in the `small/` directory to see the simulation on a smaller chunk of the landscape

- run `Rscript make_maps.R` to make the maps that SLiM uses

- `.tif` files: Yi-Ming says (3 June 2021) " The file names contain the name of the time period and the chunk of time. For example, BA13_15 means BA period that is from 15 kyr to 13 kyr."

- To run a simulation copy the required map files into a subdirectory,
    make a `params.json` JSON file with parameters that change the default,
    and run `slim ../nebria.slim` from within that directory.

## Outline of workflow:

- In each directory there is a `params.json` and a symlink to `geo_layers/`, containing the maps

- then to run a simulation with those parameter values, we:

    1. cd into that directory
    2. run `slim <relative path to nebria.slim>` in that directory, which then looks in the current
        directory for config (ie `params.json`)

- this outputs, in that directory:

    1. `sim_<seed>.trees`: the trees file (note it does this every so often, to the same file, not just at the end)
        Note that if we didn't set the seed on the slim command line (like `slim -s 123 nebria.slim`)
        then the `seed` will be almost certainly unique.
    2. `sim_<seed>.log`: a logfile of statistics
    3. `sim_<seed>_<years ago>.png`: snapshots of the simulation

- then we compute stats using `python3 compute_stats.py [name of trees file]`. This creates a subdirectory,
    which is named `<name of tree file without .trees>_stats`, and in this:

    1. The `repname` will be `{basedir}/stats_{recap_seed}_{mut_seed}_{match_rep}`,
    2. Dumps the recapitation parameters to `{repname}.json`
    3. Writes the single-location stats to `{repname}.stats.csv`
    4. Writes the pairwise statistics to `{repname}.pairstats.csv`


