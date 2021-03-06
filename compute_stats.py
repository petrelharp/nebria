import sys, os
import warnings
import time
import pyslim, tskit, msprime
import numpy as np
import pandas as pd
import json

import matplotlib.pyplot as plt
import matplotlib.collections as mc

import recap

if len(sys.argv) != 2:
    print(f"""
    Usage:
        python {sys.argv[0]} <input>.trees
    """)
    sys.exit()

ts_file = sys.argv[1]
if not os.path.isfile(ts_file):
    raise ValueError(f"File {ts_file} does not exist.")

basename = ts_file.replace(".trees", "")
basedir = f"{basename}_stats"
if not os.path.exists(basedir):
    os.makedirs(basedir)


start_time = time.time()
rng = np.random.default_rng()

# maximum distance to match observed patches with sample locations
max_dist = 10  # km
# minimum size of a simulated group to be called a 'patch'
min_patch_size = 10  # individuals
# radius within which to merge groups of individuals
patch_radius = 0.5  # km
# mutation rate
mut_rate = 2.8e-9

### number of replicates
# the total number of replicates will the product of these
replicates = {
        "recapitation": 5,
        "mutation": 5,
        "match_patch": 5,
}


print(
        f"Reading in from {ts_file}, outputting to {basename}*\n"
        f"  Merging patches of minimum size {min_patch_size} closer than {patch_radius} km of each other.\n"
        f"  Maximum distance to match patches to observed locations: {max_dist} km."
)

# SLiM doesn't write out time units yet
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# Load tree sequence, reduce to samples today 
# and remove extra population (if present)
orig_ts = pyslim.load(ts_file)
_alive_nodes = np.array(
        [n for i in orig_ts.individuals_alive_at(0)
            for n in orig_ts.individual(i).nodes
            if orig_ts.node(n).population == 1
])
orig_ts = orig_ts.simplify(_alive_nodes, keep_input_roots=True)

# then reassign these to north/middle/south
# "populations" for recapitation.

setup_demog, _ = recap.get_demography()
orig_ts = pyslim.SlimTreeSequence(recap.setup(orig_ts, setup_demog))

# keep locations to make sure these don't change through recapitation, etc
ind_locs = np.array([ind.location for ind in orig_ts.individuals()])

# Determine "patches" by merging locations within patch_radius of each other
patches = recap.assign_patches(
        orig_ts,
        patch_radius
)
print(f"Merged {orig_ts.num_individuals} into {len(patches)} patches,")

real_locs = pd.read_csv("sample_locs.csv") .rename(
        columns=lambda x: x.replace(".", "_"),
    )
real_locs.set_index('site_name', inplace=True)
real_locs.index.names = ['site_name']

# Now only keep patches of size at least min_patch_size
patches = {a:b for a, b in patches.items() if len(b) >= min_patch_size}
patch_xy = np.array([xy for xy in patches])
patch_sizes = np.array([len(patches[tuple(a)]) for a in patch_xy])

print(f"... but keeping only the {len(patches)} patches of size "
      f"at least {min_patch_size}, having a total of {sum(patch_sizes)} "
      f"individuals between them.")

# For each real location find the samples within max_dist of it,
# and the nearest patch with at least the required number of samples
real_locs["num_nearby"] = None
close_patches = { }
for name in real_locs.index:
    pdist = recap.dist(
            patch_xy,
            np.array(real_locs.loc[name, ["slim_x", "slim_y"]])
    )
    nearby = (pdist <= max_dist)
    real_locs.loc[name, "num_nearby"] = np.sum(nearby)
    close_patches[name] = np.where(nearby)[0]

####
# plot locations of individuals and patches
# and where is their nearest, large enough sampling location(s)

if True:
    # Plot all individuals
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
    sample_indivs = orig_ts.individuals_alive_at(0)
    sample_locs = np.array([orig_ts.individual(i).location[:2] for i in sample_indivs])
    ax1.scatter(sample_locs[:,0], sample_locs[:,1], marker=".", label="individuals", c='black', s=0.1)
    for ax in (ax1, ax2):
        ax.scatter(patch_xy[:,0], patch_xy[:,1], s=patch_sizes, marker="o", label="simulated samples")
        ax.scatter(real_locs["slim_x"], real_locs["slim_y"], label="real samples")
        ax.set_aspect('equal')
        ax.set_xlabel("eastings")
        ax.set_ylabel("northings")

    ax1.legend()
    plt.tight_layout()
    plt.savefig(f"{basename}.locs.png")

if True:
    # Plot patches and their matches
    fig, ax = plt.subplots(1, 1, figsize=(5, 6))
    ax.set_aspect('equal')
    ax.set_xlabel("eastings")
    ax.set_ylabel("northings")
    ax.scatter(patch_xy[:,0], patch_xy[:,1], s=patch_sizes, marker="o", label="simulated samples")
    ax.scatter(real_locs["slim_x"], real_locs["slim_y"], label="real samples", zorder=5)
    lines = []
    for k in range(real_locs.shape[0]):
        name = real_locs.index[k]
        xy0 = (real_locs["slim_x"][k], real_locs["slim_y"][k])
        for j in close_patches[name]:
            xy1 = patch_xy[j]
            lines.append([xy0, xy1])

    lc = mc.LineCollection(
            lines,
            label='matches', 
            colors=["black" if n > 0 else "red" for n in real_locs["num_nearby"]],
            linewidths=1,
            zorder=3
        )
    ax.add_collection(lc)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{basename}.locs.pdf")

if True:
    # heterozygosity at all simulated locations
    het = orig_ts.diversity([patches[tuple(xy)] for xy in patch_xy], mode='branch')

    fig, ax = plt.subplots(1, 1, figsize=(5, 6))
    ax.set_aspect('equal')
    ax.set_xlabel("eastings")
    ax.set_ylabel("northings")

    ax.scatter(patch_xy[:,0], patch_xy[:,1], s=het/10, marker="o",
            facecolors='none', edgecolors='blue',
            label="expected heterozygosity")
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"{basename}.pi.pdf")



#######
setup_time = time.time()
print(f"time: Setup done in {setup_time - start_time}. Beginning recapitation.")

for recap_rep in range(replicates["recapitation"]):

    recap_seed = rng.integers(1000000)
    recap_row = rng.integers(recap.posterior_length)
    demog, recap_params = recap.get_demography(recap_row)
    ts = msprime.sim_ancestry(
            initial_state=orig_ts,
            demography=demog,
            recombination_rate=2.48e-8,
            random_seed=recap_seed
    )

    # make sure there's no re-ordering
    for ind, loc in zip(ts.individuals(), ind_locs):
        assert np.all(loc == ind.location)

    #######
    recap_time = time.time()
    print(f"time: Recap #{recap_rep} done in {recap_time - setup_time}. Adding mutations.")

    for mut_rep in range(replicates["mutation"]):

        mut_seed = rng.integers(1000000)
        ts = msprime.sim_mutations(
                ts,
                rate=mut_rate,
                model=msprime.SLiMMutationModel(type=0),
                random_seed=mut_seed,
                keep=False,
        )

        #######
        mut_time = time.time()
        print(f"time: Mutations #{mut_rep} added in {mut_time - recap_time}. Setting up for statistics.")

        for match_rep in range(replicates["match_patch"]):

            repname = f"{basedir}/stats_{recap_seed}_{mut_seed}_{match_rep}"
            with open(f"{repname}.json", "w") as f:
                json.dump(recap_params, f)

            match_patch_vec = np.repeat(-1, len(real_locs.index))
            for k, name in enumerate(real_locs.index):
                mp = rng.choice(close_patches[name])
                ntries = 0
                while mp in match_patch_vec:
                    mp = rng.choice(close_patches[name])
                    ntries += 1
                    if ntries > 100:
                        raise ValueError("No available nearby patches to match.")
                match_patch_vec[k] = mp

            assert np.all(match_patch_vec >= 0)
            real_locs["match_patch"] = match_patch_vec

            #######
            stat_setup_time = time.time()
            print(f"time: Statistics setup in {stat_setup_time - mut_time}. Computing statistics.")


            ####
            # compute statistics

            # get matching sample sets
            sample_sets = {}
            real_locs["num_sim_samples"] = 0
            for k, r in enumerate(real_locs.itertuples()):
                row = r._asdict()
                xy = tuple(patch_xy[row["match_patch"]])
                nodes = patches[xy]
                n = row["sample_size"]
                real_locs.loc[real_locs.index[k], "num_sim_samples"] = len(nodes)
                sample_sets[row['Index']] = nodes


            # expected heterozygosity
            het_vec = np.repeat(np.nan, len(real_locs.index))
            for k in range(real_locs.shape[0]):
                name = real_locs.index[k]
                samples = sample_sets[name]
                if len(samples) > 0:
                    het_vec[k] = ts.diversity(samples, mode='site')
            real_locs["het"] = het_vec

            # write out text file
            real_locs.to_csv(f"{repname}.stats.csv")

            ###
            # pairwise stats

            pairs = pd.DataFrame(
                        np.array(
                            [[a, b] for a in real_locs.index for b in real_locs.index if a <= b]
                        ),
                        columns=["loc1", "loc2"],
            )
            pairs["dxy"] = np.nan
            pairs["Fst"] = np.nan

            has_samples = np.logical_and(
                    real_locs.loc[pairs["loc1"], "num_sim_samples"].values > 0,
                    real_locs.loc[pairs["loc2"], "num_sim_samples"].values > 0
            )

            slist = []
            nlist = []
            for n in sample_sets:
                s = sample_sets[n]
                if len(s) > 0:
                    nlist.append(n)
                    slist.append(s)

            plist = []
            for a, b in zip(pairs["loc1"][has_samples], pairs["loc2"][has_samples]):
                plist.append((nlist.index(a), nlist.index(b)))

            pairs.loc[has_samples, "dxy"] = ts.divergence(
                    sample_sets = slist,
                    indexes = plist,
                    mode = 'site',
            )
            pairs.loc[has_samples, "Fst"] = ts.Fst(
                    sample_sets = slist,
                    indexes = plist,
                    mode = 'site',
            )

            # write out text file
            pairs.to_csv(f"{repname}.pairstats.csv")

#######
end_time = time.time()
print(f"time: Total time: {end_time - start_time}.")
