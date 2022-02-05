import numpy as np
import tskit, msprime, pyslim

def dist(x, y):
    return np.sqrt(
            (x[..., 0].astype("float") - y[..., 0].astype("float")) ** 2
            + (x[..., 1].astype("float") - y[..., 1].astype("float")) ** 2
    )

def get_demography(
        conness=1000,
        selden=1000,
        army=1000,
        na=1000,
        split_time=596005.1,
):
    '''
    This is from git@github.com:yimingweng/N_ingens_ABC/scripts/model1.py
    '''
    demography = msprime.Demography()
    demography.add_population(name="SLiM", initial_size=1)
    demography.add_population(name="north", initial_size=conness) # Conness
    demography.add_population(name="center", initial_size=selden) # Selden
    demography.add_population(name="south", initial_size=army)    # Army
    demography.add_population(name="COS", initial_size=na)
    demography.add_population_split(time=split_time, derived=["north"], ancestral="COS")
    demography.add_population_split(time=split_time, derived=["center"], ancestral="COS")
    demography.add_population_split(time=split_time, derived=["south"], ancestral="COS")
    
    return demography


def setup(ts, demography):
    """
    Reassign population labels in a tree sequence so that it's ready to be
    recapitated with the three-population north/center/south model.
    This happens by assigning each node to "north" / "center" / "south"
    according to which of Conness, Selden, or Army they are closer to.
    """
    slim_coords = { # from sample_locs.csv
            "north" : np.array((62.3801745906516, 215.198439653732)), # Conness
            "center" : np.array((98.4717731879151, 142.281862598546)), # Selden
            "south" : np.array((156.664692942048, 55.1959322501973)), # Army
    }

    if ts.num_populations > 1:
        assert len(ts.samples(population=0)) == 0 # population 0 is unused
        if ts.num_populations == 3:
            # debug runs have one individual in fake population 2
            assert len(ts.samples(population=2)) == 2
        ts = ts.simplify(ts.samples(population=1), keep_input_roots=True)

    t = ts.dump_tables()
    pop_ids = {}
    pop_names = [p.name for p in demography.populations]
    for x in slim_coords:
        assert x in pop_names, f"'{x}' is not in the demography's population names"
    t.populations.clear()
    for j, pn in enumerate(pop_names):
        if pn == "SLiM":
            # keep metadata, including SLiM id for existing population
            pop_md = ts.population(0).metadata
        else:
            pop_md = pyslim.default_slim_metadata("population")
            pop_md.update({ "slim_id": j + 2 })
        pop_md.update({ "name": pn })
        k = t.populations.add_row(
                metadata=pop_md
        )
        pop_ids[pn] = k
    
    ind_locs = t.individuals.location.reshape((t.individuals.num_rows, 3))
    dist_pop_names = ["north", "center", "south"]
    dist_locs = np.column_stack([
            dist(slim_coords[k], ind_locs)
            for k in dist_pop_names
    ])
    remap = np.array([pop_ids[pn] for pn in dist_pop_names])
    # which_pop[i] is the id of the population to which individual i is the closest
    which_pop = remap[
            np.argmin(dist_locs, axis=1) 
    ]

    max_t = max(t.nodes.time)
    t.nodes.clear()
    for n in ts.nodes():
        new_pop = pop_ids["SLiM"]
        if n.time == max_t and n.individual != tskit.NULL:
            new_pop = which_pop[n.individual]
        t.nodes.append(n.replace(population=new_pop))

    ts = t.tree_sequence()
    return ts


