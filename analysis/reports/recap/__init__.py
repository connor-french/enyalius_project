import numpy as np
import tskit, msprime, pyslim
import importlib.resources
from . import data

# this is the way to include a data file, I guess?
# https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package
_post_file = importlib.resources.open_text(data, 'posterior.txt')
_header = _post_file.readline().strip().split()
_post = np.loadtxt(_post_file, skiprows=1)
# keys are ['T2', 'T1', 'CS', 'AS', 'NE', 'Na', 'Nc', 'Ns']
posterior = { _header[k]: _post[:,k] for k in range(len(_header)) }
posterior_length = _post.shape[0]


def dist(x, y):
    return np.sqrt(
            (x[..., 0].astype("float") - y[..., 0].astype("float")) ** 2
            + (x[..., 1].astype("float") - y[..., 1].astype("float")) ** 2
    )


def get_demography(rep=0):
    '''
    This is from git@github.com:yimingweng/N_ingens_ABC/scripts/model2.py
    '''
    params = {
        k: posterior[k][rep] for k in posterior
    }

    demography = msprime.Demography()
    demography.add_population(name="SLiM", initial_size=1)
    demography.add_population(name="north", initial_size=params["Nc"]) # Conness
    demography.add_population(name="center", initial_size=params["Ns"]) # Selden
    demography.add_population(name="south", initial_size=params["Na"])    # Army
    demography.add_population(name="COS", initial_size=params["NE"])
    demography.add_admixture(
            time=params["T2"],
            derived="center",
            ancestral=["north", "south"],
            proportions=[params["CS"], 1 - params["CS"]]
    )
    demography.add_population_split(time=params["T1"], derived=["north", "south"], ancestral="COS")
    return demography, params


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


def assign_patches(ts, patch_radius):
    """
    Merge individuals into "patches" using a simple greedy algorithm.
    The output is a dictionary that is indexed by (x, y) tuples of spatial
    location, and values are lists of the individual IDs assigned to that patch.

    Strategy:
    1. Round locations to a resolution smaller than `patch_radius`.
    2. Sort unique locations by number of individuals.
    3. Merge any location with the largest other location within `patch_radius`
       that is not itself merged to a larger location.
    """
    sample_indivs = ts.individuals_alive_at(0)
    sample_locs = np.array([ts.individual(i).location[:2] for i in sample_indivs])

    # first round to this many binary digits
    digits2 = 2 + int( (-1) * np.floor(np.log2(patch_radius)) )
    def coarsen(x):
        return np.round( (2 ** digits2) * x ) / (2 ** digits2)

    patches = {}
    for i, (x, y) in enumerate(sample_locs):
        xy = (coarsen(x), coarsen(y))
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
            pdist = dist(xy, _sorted_patches[(j+1):,:2])
            for k, d in zip(range(j+1, len(_sorted_patches)), pdist):
                xy1 = (_sorted_patches[k, 0], _sorted_patches[k, 1])
                if ((d <= patch_radius) and
                        (xy1 not in _merge_with)):
                    _merge_with[xy1] = tuple(xy)

    for xy1, xy0 in _merge_with.items():
        patches[xy0].extend(patches[xy1])
        del patches[xy1]

    return patches
