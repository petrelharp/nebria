import sys, os
import warnings
import time
import pyslim, tskit, msprime
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.collections as mc

import recap

start_time = time.time()

rng = np.random.default_rng()

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

# maximum distance to match observed patches with sample locations
max_dist = 10  # km
# minimum size of a simulated group to be called a 'patch'
min_patch_size = 10  # individuals
# radius within which to merge groups of individuals
patch_radius = 0.4  # km

# mutation rate
mut_rate = 2.8e-9

print(
        f"Reading in from {ts_file}, outputting to {basename}*\n"
        f"  Merging patches of minimum size {min_patch_size} closer than {patch_radius} km of each other.\n"
        f"  Maximum distance to match patches to observed locations: {max_dist} km."
)

# SLiM doesn't write out time units yet
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# recapitate and mutate
ts = tskit.load(ts_file)

# Reduce to samples today, then reassign these to north/middle/south
# "populations" for recapitation.
ts = ts.simplify(ts.samples(population=1, time=0), keep_input_roots=True)

demog = recap.get_demography()
ts = recap.setup(ts, demog)

#######
setup_time = time.time()
print(f"time: Setup done in {setup_time - start_time}. Beginning recapitation.")

ts = msprime.sim_ancestry(
        initial_state=ts,
        demography=demog,
        recombination_rate=2.48e-8,
)

#######
recap_time = time.time()
print(f"time: Recap done in {recap_time - setup_time}. Adding mutations.")

ts = msprime.sim_mutations(
        ts,
        rate=mut_rate,
        model=msprime.SLiMMutationModel(type=0)
)

#######
mut_time = time.time()
print(f"time: Mutations added in {mut_time - recap_time}. Setting up for statistics.")

ts = pyslim.SlimTreeSequence(ts)

real_locs = pd.read_csv("sample_locs.csv") .rename(
        columns=lambda x: x.replace(".", "_"),
    )
real_locs.set_index('site_name', inplace=True)
real_locs.index.names = ['site_name']

# Merge locations within patch_radius of each other
patches = recap.assign_patches(ts, patch_radius)
print(f"Merged {ts.num_individuals} into {len(patches)} patches,")

# Now only keep patches of size at least min_patch_size
patches = {a:b for a, b in patches.items() if len(b) >= min_patch_size}
patch_xy = np.array([xy for xy in patches])
patch_sizes = np.array([len(patches[tuple(a)]) for a in patch_xy])

print(f"  ... but keeping only the {len(patches)} patches of size "
      f"at least {min_patch_size}, having a total of {sum(patch_sizes)} "
      f"individuals between them.")

# For each real location find the samples within max_dist of it,
# and the nearest patch with at least the required number of samples
real_locs["num_nearby"] = None
real_locs["closest_patch"] = -1  # BEWARE!!! but pandas has no reasonable missing data type
close_patches = { }
for k in range(real_locs.shape[0]):
    name = real_locs.index[k]
    pdist = recap.dist(
            patch_xy,
            np.array(real_locs.loc[real_locs.index[k], ["slim_x", "slim_y"]])
    )
    nearby = (pdist <= max_dist)
    real_locs.loc[real_locs.index[k], "num_nearby"] = np.sum(nearby)
    close_patches[name] = np.where(nearby)[0]
    # only match to patches with at least as many samples as we have in the real data
    big_enough = np.where(patch_sizes >= real_locs["sample_size"][k])[0]
    real_locs.loc[real_locs.index[k], "closest_patch"] = big_enough[np.argmin(pdist[big_enough])]
    assert pdist[real_locs.loc[real_locs.index[k], "closest_patch"]] <= max_dist 


####
# plot locations of individuals and patches
# and where is their nearest, large enough sampling location(s)

# all individuals
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
sample_indivs = ts.individuals_alive_at(0)
sample_locs = np.array([ts.individual(i).location[:2] for i in sample_indivs])
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

# patches and their matches
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

#######
stat_setup_time = time.time()
print(f"time: Statistics setup in {stat_setup_time - mut_time}. Computing statistics.")


####
# compute statistics


# get matching sample sets
sample_sets = {}
_sampled = set()
real_locs["num_sim_samples"] = 0
for k, r in enumerate(real_locs.itertuples()):
    row = r._asdict()
    xy = tuple(patch_xy[row["closest_patch"]])
    nodes = list(set(patches[xy]) - set(_sampled))
    n = row["sample_size"]
    if len(nodes) < n:
        print(f"Not enough samples for {row['Index']} left: only {len(nodes)} < {n}.")
    s = rng.choice(nodes, size=min(len(nodes), n), replace=False)
    real_locs.loc[real_locs.index[k], "num_sim_samples"] = len(s)
    sample_sets[row['Index']] = s
    _sampled.update(s)


# expected heterozygosity
real_locs["het"] = np.nan
for k in range(real_locs.shape[0]):
    name = real_locs.index[k]
    samples = sample_sets[name]
    if len(samples) > 0:
        pi = ts.diversity(samples, mode='site')
        real_locs.loc[real_locs.index[k], 'het'] = pi


# write out text file
real_locs.to_csv(f"{basename}.stats.csv")


# heterozygosity at all simulated locations
het = ts.diversity([patches[tuple(xy)] for xy in patch_xy], mode='branch')

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
pairs.to_csv(f"{basename}.pairstats.csv")

#######
end_time = time.time()
print(f"time: Statistics computed in {end_time - stat_setup_time}.")
print(f"time: Total time: {end_time - start_time}.")
