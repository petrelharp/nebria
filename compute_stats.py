import sys, os
import warnings
import pyslim, tskit, msprime
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.collections as mc

import recap

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

# maximum distance to match samples out to
max_dist = 10  # km
# minimum size of a simulated group to be called a 'patch'
min_patch_size = 10  # individuals
# radius within which to merge groups of individuals
patch_radius = 0.05  # km

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

ts = msprime.sim_ancestry(
        initial_state=ts,
        demography=demog,
        recombination_rate=2.48e-8,
)

ts = msprime.sim_mutations(
        ts,
        rate=mut_rate,
        model=msprime.SLiMMutationModel(type=0)
)

ts = pyslim.SlimTreeSequence(ts)

real_locs = pd.read_csv("sample_locs.csv") .rename(
        columns=lambda x: x.replace(".", "_"),
    )
real_locs.set_index('site_name', inplace=True)
real_locs.index.names = ['site_name']

sample_indivs = ts.individuals_alive_at(0)
sample_locs = np.array([ts.individual(i).location[:2] for i in sample_indivs])

##########
# Find the distinct patches of simulated individuals:
# `patches[xy]` will contain a list of the individual IDs
#    at the patch with coordinates xy
#
# Strategy:
# 1. Round locations to the nearest 10m
# 2. Sort unique locations by number of individuals
# 3. Merge any location with the largest other location within max_distance
#    that is not itself merged to a larger location.
# 4. Retain merged locations with at least min_patch_size individuals.


patches = {}
for i, (x, y) in enumerate(sample_locs):
    xy = (round(x, 1), round(y, 1))
    if xy not in patches:
        patches[xy] = []
    patches[xy].append(i)

# now merge very close ones:
_sorted_patches = list(patches.keys())
# these will be in decreasing order of size
_sorted_patches.sort(key = lambda x: len(patches[x]), reverse=True)
_sorted_patches = np.array(_sorted_patches)
_merge_with = {}
for j, xyn in enumerate(_sorted_patches):
    xy = xyn[:2]
    if tuple(xy) not in _merge_with:
        pdist = recap.dist(xy, _sorted_patches[(j+1):,:2])
        for k, d in zip(range(j+1, len(_sorted_patches)), pdist):
            xy1 = (_sorted_patches[k, 0], _sorted_patches[k, 1])
            if ((d <= patch_radius) and
                    (xy1 not in _merge_with)):
                _merge_with[xy1] = tuple(xy)

for xy1, xy0 in _merge_with.items():
    patches[xy0].extend(patches[xy1])
    del patches[xy1]

# now only keep patches of size at least min_patch_size
patches = {a:b for a, b in patches.items() if len(b) >= min_patch_size}
patch_xy = np.array([xy for xy in patches])
patch_sizes = np.array([len(patches[tuple(a)]) for a in patch_xy])

# for each real location find the samples within max_dist of it
real_locs["num_nearby"] = None
real_locs["close_patch"] = -1  # BEWARE!!! but pandas has no reasonable missing data type
for k in range(real_locs.shape[0]):
    pdist = recap.dist(
            patch_xy,
            np.array(real_locs.loc[real_locs.index[k], ["slim_x", "slim_y"]])
    )
    nearby = np.where(pdist <= max_dist)[0]
    real_locs.loc[real_locs.index[k], "num_nearby"] = len(nearby)
    real_locs.loc[real_locs.index[k], "close_patch"] = np.argmin(pdist)


####
# plot locations of individuals
# and where is their nearest, largest sampling location(s)

fig, ax = plt.subplots(1, 1, figsize=(5, 6))
ax.set_aspect('equal')
ax.set_xlabel("eastings")
ax.set_ylabel("northings")

lc = mc.LineCollection(
        [[(real_locs["slim_x"][k], real_locs["slim_y"][k]),
            patch_xy[real_locs["close_patch"][k]]]
            for k in np.where(real_locs["close_patch"] >= 0)[0]],
        label='matches', 
        colors=["black" if n > 0 else "red" for n in real_locs["num_nearby"]],
        linewidths=1,
        zorder=0
    )
ax.add_collection(lc)
ax.scatter(patch_xy[:,0], patch_xy[:,1], s=patch_sizes, marker="o", label="simulated samples")
ax.scatter(real_locs["slim_x"], real_locs["slim_y"], label="real samples")

ax.legend()

plt.tight_layout()
plt.savefig(f"{basename}.locs.pdf")


####
# compute statistics


# get matching sample sets
sample_sets = {}
_sampled = set()
real_locs["num_sim_samples"] = 0
for k, r in enumerate(real_locs.itertuples()):
    row = r._asdict()
    xy = tuple(patch_xy[row["close_patch"]])
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
