import pyslim, tskit, msprime
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.collections as mc

rng = np.random.default_rng()

ts_file = "nebria_3541435757397.trees"
basename = ts_file.replace(".trees", "")

# maximum distance to match samples out to
max_dist = 10  # km
# minimum size of a simulated group to be called a 'patch'
min_patch_size = 10  # individuals
# radius within which to merge groups of individuals
patch_radius = 0.05  # km

# mutation rate
mut_rate = 1e-8

# recapitate and mutate
ts = pyslim.load(ts_file)
ts = pyslim.recapitate(ts, ancestral_Ne=1e4, recombination_rate=1e-8)
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

# subpopulation 2, if present, is our fake population for visualization
sample_indivs = []
for i in ts.individuals_alive_at(0):
    ind = ts.individual(i)
    ns = [ts.node(n).is_sample() for n in ind.nodes]
    if ((ind.metadata['subpopulation'] == 1)
            and (sum(ns) == len(ns))):
        sample_indivs.append(i)

sample_indivs = np.array(sample_indivs)
sample_locs = np.array([ts.individual(i).location[:2] for i in sample_indivs])

def dist(x, y):
    return np.sqrt(
            (x[..., 0].astype("float") - y[..., 0].astype("float")) ** 2
            + (x[..., 1].astype("float") - y[..., 1].astype("float")) ** 2
    )

##########
# Find the distinct patches of simulated individuals:
# `patches[xy]` will contain a list of the individual IDs
#    at the patch with coordinates xy

patches = {}
for i, (x, y) in enumerate(sample_locs):
    j = (round(x, 2), round(y, 2))
    if j not in patches:
        patches[j] = []
    patches[j].append(i)

# now merge very close ones:
_sorted_patches = list(patches.keys())
_sorted_patches.sort(key = lambda x: len(patches[x]), reverse=True)
_sorted_patches = np.array(_sorted_patches)
_merge_with = {}
for j, xy in enumerate(_sorted_patches):
    if tuple(xy) not in _merge_with:
        pdist = dist(xy, _sorted_patches[(j+1):,:2])
        for k, d in zip(range(j+1, len(_sorted_patches)), pdist):
            xy1 = (_sorted_patches[k, 0], _sorted_patches[k, 1])
            if ((d <= patch_radius) and
                    (xy1 not in _merge_with)):
                _merge_with[xy1] = tuple(xy)

for xy1, xy0 in _merge_with.items():
    patches[xy0].extend(patches[xy1])
    del patches[xy1]

# now only keep patches of size at least a certain size
patches = {a:b for a, b in patches.items() if len(b) >= min_patch_size}
patch_xy = np.array([xy for xy in patches])
patch_sizes = np.array([len(patches[tuple(a)]) for a in patch_xy])

# for each real location find the samples within max_dist of it
real_locs["num_nearby"] = None
real_locs["close_patch"] = -1  # BEWARE!!! but pandas has no reasonable missing data type
for k in range(real_locs.shape[0]):
    pdist = dist(
            patch_xy[:,:2],
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
