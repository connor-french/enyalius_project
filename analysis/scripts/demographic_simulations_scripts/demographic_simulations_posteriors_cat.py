###### Exploratory simulations to see impacts of parameter values on summary statistics #######

# import packages
from os.path import isfile
import numpy as np
from numpy import ma
import random
import msprime
import rasterio
from rasterio.mask import mask
import allel
import scipy.stats as stats
import pandas as pd
import itertools
from math import radians
from sklearn.metrics.pairwise import haversine_distances
from sklearn import linear_model
from shapely.geometry import Point
from esda.moran import Moran
from libpysal.weights import Voronoi
import time
from scipy.interpolate import NearestNDInterpolator
import numba
import tracemalloc
from scipy.stats import entropy
import warnings
from multiprocessing import Pool
from functools import partial


# ignoring the performance warnings from how I build the pandas dataframes
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# initialize MPI
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()

#### Define functions #####

# return nan if it's trying to take the mean or std of something with no values
def nanmean(lst):
    try:
        if len(lst) == 0:
            return np.nan
        return np.nanmean(lst)
    except:
        return np.nan
    
def nanstd(lst):
    try:
        if len(lst) == 0:
            return np.nan
        return np.nanstd(lst)
    except:
        return np.nan

## raster_to_k ##
def raster_to_k(raster, transformation = "linear", max_local_k = 100, normalize = False, threshold = None, inflection_point = 0.5, slope = 0.05):

    d = raster.filled(0)

    if normalize:
        def normalize(rast):
            return (rast - np.min(rast)) / (np.max(rast) - np.min(rast))
        
        d = normalize(d)
    
    if transformation == "linear":
        t = d * max_local_k
        t = np.ceil(t)

    if transformation == "hinge":
        t = d
        t[t >= threshold] = max_local_k
        t[t < threshold] = 1e-10

    if transformation == "sigmoid":
        def sigmoid(x, a, b):
            y = 1.0 / (1.0 + np.ma.exp(-(x - a) / b))
            return(y)
        
        sigmoid_v = np.vectorize(sigmoid)
        t = sigmoid_v(d, inflection_point, slope) * max_local_k
        t = np.ceil(t)
    
    t[t == 0] = 1e-10
    
    return(t)
  
## migration_matrix ##
def migration_matrix(populations, rate, scale=True):
    d = populations.shape[0] * populations.shape[1]
    M = np.zeros((d, d))

    for i in range(populations.shape[0]):
        for j in range(populations.shape[1]):
            current_index = i * populations.shape[1] + j
            # check the neighboring populations and calculate the migration rates. Looping through tuples was a neat trick! Saved a lot of time.
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                # check for edges
                if 0 <= i + di < populations.shape[0] and 0 <= j + dj < populations.shape[1]:
                    neighbor_index = (i + di) * populations.shape[1] + (j + dj)
                    if scale:
                        # mig = donor / recipient * rate unless the pop size is zero
                        M[current_index, neighbor_index] = (populations[i , j] / populations[i + di, j + dj]) * rate if populations[i + di, j + dj] > 1e-9 and populations[i, j] > 1e-9 else 0
                    else:
                        M[current_index, neighbor_index] = rate if populations[i + di, j + dj] > 1e-9 and populations[i, j] > 1e-9 else 0


    return M

## stepping_stone2d ##
def stepping_stone2d(deme_size_matrix, rate, scale=True):
    assert len(deme_size_matrix.shape) <= 3

    n, m = deme_size_matrix.shape
    N = n * m
    model = msprime.Demography.isolated_model(deme_size_matrix.reshape(N))

    # set population names
    for j in range(n):
        for k in range(m):
            index = j * m + k
            model.populations[index].name = f"pop_{j}_{k}"
    
    # setup migration rate matrices
    if np.array(rate).ndim == 0:
        if scale:
            model.migration_matrix = migration_matrix(deme_size_matrix, rate, scale=True)
        else: 
            model.migration_matrix = migration_matrix(deme_size_matrix, rate, scale=False)
    else:
        assert rate.shape == (N, N), f"Expected a migration matrix with the shape {(N, N)} and instead got {rate.shape}"
        model.migration_matrix = rate

    
    return model

## add_landscape_change ##
def add_landscape_change(model, k_stack, timestep = 1, rate = 0.001, scale=True):
    # iterate through the first dimension of a 3D array, where the array represents different time steps of population size change
    # omit the most ancient time step (have to set its migration rate differently)
    for step in range(1, k_stack.shape[0] - 1):
        # get the population size values of the current array
        kmat = k_stack[step]

        # get the population size values of array from the more ancient time step
        kmat_anc = k_stack[step + 1]

        # get the population size values array from the more recent time step
        kmat_prev = k_stack[step - 1]

        # get the shape of the array
        n, m = kmat.shape

        ##### Update population sizes #####
        # add population size changes according to the values of the current array
        for j in range(n):
            for k in range(m):
                # only update the population size if it is different from the previous time point
                if kmat[j, k] != kmat_prev[j, k]:
                    # add a demographic change to each cell in the raster
                    model.add_population_parameters_change(time=step * timestep, population=f"pop_{j}_{k}", initial_size=kmat[j, k])

        ##### Update migration rates #####
        # add migration rate change for each time step
        # this is updating time steps from the present to the past
        # ## iterate through the population sizes
        for i in range(n):
            for j in range(m):
                ## also index the neighboring cells
                for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    ## check for edges
                    if 0 <= i + di < kmat.shape[0] and 0 <= j + dj < kmat.shape[1]:
                        ## only update migration if the donor and recipient population sizes are different between time points
                        if kmat_prev[i + di, j + dj] != kmat[i + di, j + dj] and kmat[i, j] != kmat_prev[i, j]:
                            if scale:
                                ## mig = donor / recipient * rate unless the pop size is zero
                                r = (kmat[i, j]  / kmat[i + di, j + dj]) * rate if kmat[i + di, j + dj] > 1e-9 and kmat[i, j] > 1e-9 else 0
                                model.add_migration_rate_change(time = step * timestep, rate = r, source=f"pop_{i}_{j}", dest=f"pop_{i + di}_{j + dj}")
                            else:
                                model.add_migration_rate_change(time = step * timestep, rate = rate, source=f"pop_{i}_{j}", dest=f"pop_{i + di}_{j + dj}")
                        elif kmat_prev[i + di, j + dj] != kmat[i + di, j + dj] and kmat[i, j] != kmat_prev[i, j] and kmat_anc[i, j] <= 1e-9:
                            ## have the deme migrate to neighbors if the more ancient time step has an empty deme
                            if scale:
                                r = (kmat[i, j] / kmat[i + di, j + dj]) * rate if kmat[i + di, j + dj] > 1e-9 and kmat[i, j] > 1e-9 else 0
                                model.add_migration_rate_change(time = step * timestep, rate = r, source=f"pop_{i}_{j}", dest=f"pop_{i + di}_{j + dj}")
                            else:
                                model.add_migration_rate_change(time = step * timestep, rate = rate, source=f"pop_{i}_{j}", dest=f"pop_{i + di}_{j + dj}")
    
    return model


# function to split a landscape according to the "STRUCTURE" population assigned to sampled individuals (these populations can be considered admixture populations)
# raster_path is the path to the landscape raster (i.e. SDM projection) that you want to divide. 
# coords is a list of tuples of the coordinates assigned to each individual in the empirical data set
# admix_id is the population ID assigned to each individual in the empirical dataset. For example, if STRUCTURE finds a "best K" of two, every individual in the data set is assigned to one of two populations. There should be one ID for every coordinate pair.
# outfile is the path to where you want to write the new raster, including the ".tif" extension to write a GeoTIFF. Right now I'm just supporting the GeoTIFF format to keep things simple. 
# band_index is the index of the raster you want to read in, in case your raster contains multiple bands. The default is one (note- rasterio begins indexing at 1 for raster bands)
# mask_rast is a boolean of whether you want to mask the interpolation by the SDM landscape or let it fill the entire extent. For the purposes of parameterizing the 2DSS, you want it to fill the entire extent, but it is useful to visualize it on the landscape. 
def split_landscape_by_admix_pop(raster_path, coords, admix_id, outfile, band_index=1, mask_rast=False):

    # open the raster
    r = rasterio.open(raster_path)
    # read in the raster. This imports the raster as a numpy array
    r2 = r.read(band_index, masked=True)

    # get the x,y indices of the empirical sampled cells
    indices_x = []
    indices_y = []
    for xy in coords:
    
        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]
        
        # mask the raster with the points
        out_image = mask(r, pt2, nodata="nan", filled = False)

        # oi first raster
        oi = out_image[0][0]

        # get the locality index
        indices_tup = ma.nonzero(oi)

        indices_x.append(indices_tup[0])
        indices_y.append(indices_tup[1])

    indices_x = np.concatenate(indices_x)
    indices_y = np.concatenate(indices_y)

    # get all of the x, y indices of the input raster
    r_x, r_y = np.indices(r.shape)

    interp = NearestNDInterpolator(list(zip(indices_x, indices_y)), admix_id)
    z = interp(r_x, r_y)

    # apply mask of the SDM landscape
    if mask_rast:
         z = ma.masked_array(z, r2.mask, fill_value=-9, dtype=float)

    # check to make sure the filepath contains the ".tif" suffix, then write the file out
    try: 
        root, ext = splitext(outfile)
        if ext != ".tif" or "~" in root:
            raise SyntaxError
        print(f"Landscape with K assignment categories written to {outfile}.")
    except SyntaxError:
        print("The outfile cannot use tilde (~) path expansion and must have a .tif extension.")

    
    with rasterio.open(
            outfile,
            mode="w",
            driver="GTiff",
            height=z.shape[0],
            width=z.shape[1],
            count=1,
            dtype=z.dtype,  
            crs=r.crs,
            transform=r.transform,
            nodata=-9
        ) as new_dataset:
            new_dataset.write(z, 1)
    

# admix_id_rast is the raster that contains which admixture populations (defined in a "STRUCTURE" type analysis) each deme belongs to. It assumes that the populations are numbered sequentially, starting at 1. 
# ancestral_size_list is a list of ancestral population sizes that correspond with each admixture population ID
# ancestral_merge_time_list is a list of times when the ancestral populations merge with each other. It should be either of size N - 1, where N is the number of ancestral populations, or 1, where all ancestral populations merge at the same time in the past. The assumed branching pattern is from left to right (e.g. for three ancestral populations, the pattern is ((ANC_1, ANC_2), ANC_3))
# ancestral_merge_size_list is a list of the sizes of the new ancestral populations after the merge. It should be the same size as the ancestral_merge_time_list
def add_nonspatial_phase(model, ancestral_size_list, merge_time, admix_id_rast=None, ancestral_merge_time_list=None, ancestral_merge_size_list=None):

    if admix_id_rast is None:
        # add an ancestral population
        model.add_population(name = "ANC", initial_size=ancestral_size_list[0])

        # get names of populations for initiating the merge
        pop_names = [[pop.name] for pop in model.populations if "ANC" not in pop.name]
        
        # add the time when the spatial simulation collapses into the collecting phase
        [model.add_population_split(time = merge_time, derived = name, ancestral = "ANC") for name in pop_names]
    else:

        # make sure the admixture ID raster has the same number of populations as the raster used in the demographic modeling
        try:
            admix_id_1d = admix_id_rast.ravel()
            if len(model.populations) != len(admix_id_1d):
                raise ValueError
        except ValueError:
            print(f"The number of demes in the demographic model is {len(model.populations)}, while the number of demes in the admixture ID raster is {len(admix_id_1d)}. They need to be the same.")

        
        # add a new ancestral population for each admixture population
        for i in range(1, len(ancestral_size_list) + 1):
            anc_pop_name = f"ANC_{i}"
            model.add_population(name = anc_pop_name, initial_size=ancestral_size_list[i - 1])

        # merge each deme in the simulation into its respective ancestral population
        for i in range(len(model.populations)):
            pop_name = model.populations[i].name
            if "ANC" not in pop_name:
                anc_pop = f"ANC_{admix_id_1d[i]}"
                model.add_population_split(time = merge_time, derived = [pop_name], ancestral = anc_pop)

        # merge ancestral populations at their respective times, if we want that behavior
        if ancestral_merge_time_list is not None:
            try:
                if len(ancestral_merge_time_list) != len(ancestral_size_list) + 1 and len(ancestral_merge_time_list) != 1:
                    raise ValueError
            except ValueError:
                print("The ancestral merge list should be of length N - 1 or 1")
            
            # make sure the ancestral merge list is either of size 1 or N - 1
            if len(ancestral_size_list) > 1 and len(ancestral_merge_time_list) > 1:
                for i in range(1, len(ancestral_size_list)):
                    if i == 1:
                        n1 = f"ANC_{i}"
                        n2 = f"ANC_{i + 1}"
                        anc_n = f"ANC_{i}_{i + 1}"
                    else:
                        anc_num_str = "_".join(map(str, range(1, i + 1)))
                        n1 = f"ANC_{anc_num_str}"
                        n2 = f"ANC_{i + 1}"
                        anc_n = f"ANC_{anc_num_str}_{i + 1}"
                    ## add the new ancestral populations to the simulation
                    model.add_population(name = anc_n, initial_size=ancestral_merge_size_list[i - 1])
                    model.add_population_split(time = ancestral_merge_time_list[i - 1], derived = [n1, n2], ancestral = anc_n)
            ## if there are multiple ancestral populations, but a single merge time
            elif len(ancestral_size_list) > 1 and len(ancestral_merge_time_list) == 1:
                # get a list of the ancestral population names
                anc_der_pops = []
                for i in range(1, len(ancestral_size_list) + 1):
                    anc_der_name = f"ANC_{i}"
                    anc_der_pops.append(anc_der_name)
                # make the name of the most ancestral population
                anc_num_str = "_".join(map(str, range(1, len(ancestral_size_list) + 1)))
                anc_n = f"ANC_{anc_num_str}"
                ## add the new ancestral populations to the simulation
                model.add_population(name = anc_n, initial_size=ancestral_merge_size_list[0])
                model.add_population_split(time = ancestral_merge_time_list[0], derived = anc_der_pops, ancestral = anc_n)


    return model

# Create a dictionary to map deme IDs to their respective admixture population assignments. Requires output sample dictionary from `coords_to_sample_dict_empirical()` and a raster with admixture population IDs for each deme, output from `split_landscape_by_admix_pop()`.  
# Returns a dictionary of the form `{admix_id: deme_id}`
def admix_to_deme_dict(admix_rast_path, pops_dict):
    
    # read in admixture raster
    with rasterio.open(admix_rast_path, "r") as src:
        r=src.read()

    # unravel the raster into one dimension for indexing
    r_1d = r.ravel()

    admix_dict = {}
    for key in pops_dict:
        k = r_1d[key]
        if k in admix_dict.keys():
            admix_dict[k].append(key)
        else:
            admix_dict[k] = [key]

    return admix_dict


## max_thresh_from_coords ##

def max_thresh_from_coords(raster, coordinates):

    xy = [Point(xy) for xy in coordinates]
    # mask the raster, so only cells that the coordinates overlap are masked
    out_image = mask(raster, xy, nodata="nan", filled = False)

    # find the minimum value of the masked raster
    max_thresh = np.min(out_image[0])

    return max_thresh


# setup demography
# returns a tuple of a Demography object and a summary of population sizes
def setup_demography(raster, transformation, max_k, admix_rast, merge_time, ancestral_size_list = None, ancestral_merge_size_list=None, ancestral_merge_time_list=None, mig_rate = 0.001, thresh=0, tstep = 1000):

    k = raster_to_k(raster=raster, transformation=transformation, threshold = thresh, max_local_k=max_k)

    m = stepping_stone2d(k[0], rate = mig_rate)

    m = add_landscape_change(model = m, k_stack = k, rate=mig_rate, timestep = tstep)

    # for my Enyalius simulations, just making the ancestral size the size of the largest admixture population
    # I don't want to change the underlying code so the more flexible option is still available
    ancestral_merge_size_list = [max(ancestral_size_list)]
    m = add_nonspatial_phase(m, 
                         ancestral_size_list=ancestral_size_list, 
                         merge_time=merge_time, 
                         admix_id_rast=admix_rast, 
                         ancestral_merge_size_list=ancestral_merge_size_list, 
                         ancestral_merge_time_list=ancestral_merge_time_list)

    m.sort_events()

    ##############################################
    ############ DEMOGRAPHY SUMMARY ##############
    ##############################################

    # DataFrame of population demography summary
    ## summarize individual numbers for each landscape in the time series
    total_inds_list = []
    deme_size_list = []
    num_occupied_list = []
    for tstep in k:
        tstep_nz = tstep[tstep > 1e-9] 
        ### total number of individuals across the landscape
        total_inds = np.sum(tstep_nz)
        ### average deme size for nonzero demes
        avg_nonzero_deme_size = nanmean(tstep_nz)
        ### number of occupied demes across the landscape
        num_occ = len(tstep_nz.ravel())

        total_inds_list.append(total_inds)
        deme_size_list.append(avg_nonzero_deme_size)
        num_occupied_list.append(num_occ)

    total_inds_list = np.array(total_inds_list)
    deme_size_list = np.array(deme_size_list)
    num_occupied_list = np.array(num_occupied_list)


    ## total number of individuals across the landscape in the present landscape
    total_inds_curr = total_inds_list[0]
    ## average number of individuals per deme in demes that contain individuals in the present landscape
    avg_deme_size_curr = deme_size_list[0]
    ## total number of individuals across the landscape in the largest landscape and the index of the largest landscape
    total_inds_largest = np.max(total_inds_list)
    total_inds_largest_index = np.argmax(total_inds_list)

    ## average number of individuals across the landscape in the landscape with the largest average size and its index
    avg_deme_size_largest = np.max(deme_size_list)
    avg_deme_size_largest_index = np.argmax(deme_size_list)

    ## average and standard deviation of total individuals across the landscape across timesteps
    total_inds_avg = nanmean(total_inds_list)
    total_inds_sd = nanstd(total_inds_list)

    ## average and standard deviation of average number of average number of individuals per deme
    avg_deme_size_avg = nanmean(deme_size_list)
    avg_deme_size_sd = nanstd(deme_size_list)

    ## number of occupied demes in the current landscape
    num_occ_curr = num_occupied_list[0]
    ## number of occupied demes in the landscape with the largest number of occupied demes and its index
    num_occ_largest = np.max(num_occupied_list)
    num_occ_largest_index = np.argmax(num_occupied_list)
    ## average number of occupied demes across timesteps
    num_occ_average = nanmean(num_occupied_list)
    ## standard deviation of the number of occupied demes
    num_occ_sd = nanstd(num_occupied_list)

    # pandas DataFrame of the demography summaries
    demo_df = pd.DataFrame(
            data = {
                "total_inds_curr": total_inds_curr,
                "avg_deme_size_curr": avg_deme_size_curr,
                "total_inds_largest": total_inds_largest,
                "total_inds_largest_index": total_inds_largest_index,
                "avg_deme_size_largest": avg_deme_size_largest,
                "avg_deme_size_largest_index": avg_deme_size_largest_index,
                "total_inds_avg": total_inds_avg,
                "total_inds_sd": total_inds_sd,
                "avg_deme_size_avg": avg_deme_size_avg,
                "avg_deme_size_sd": avg_deme_size_sd,
                "num_occ_curr": num_occ_curr,
                "num_occ_largest": num_occ_largest,
                "num_occ_largest_index": num_occ_largest_index,
                "num_occ_average": num_occ_average,
                "num_occ_sd": num_occ_sd
            }, index = [0]
        )
    
    for i in range(len(ancestral_size_list)):
        demo_df[f"admix_n_a{i + 1}"] = ancestral_size_list[i]
    
    for i in range(len(ancestral_merge_size_list)):
        demo_df[f"admix_anc_n_a{i + 1}"] = ancestral_merge_size_list[i]
    
    for i in range(len(ancestral_merge_time_list)):
        demo_df[f"anc_merge_time_a{i + 1}"] = ancestral_merge_time_list[i]

    demo_df = np.round(demo_df, 6)

    return m, demo_df


# m is the demography
def sim_ts(m, recomb_rate, seqlen = 1e6, pops_dict=None, multi_merge = False, duration = None, randseed=1):
# if there is no user-supplied pops_dict, sample all populations with at least two individuals
    if pops_dict is None:
        # get populations with at least two individuals to sample for coalescent simulations
        nzp = []
        
        for i in range(len(m.populations)):
            if m.populations[i].initial_size >= 2.0 and m.populations[i].name != "ANC":
                nzp.append(i)

        pops_dict = {key: 2 for key in nzp}

    if multi_merge:
        ts = msprime.sim_ancestry(pops_dict, 
                                sequence_length=seqlen, 
                                demography = m, 
                                model=[
                                    msprime.DiracCoalescent(duration=duration, psi = 0.1, c = 1),
                                    msprime.StandardCoalescent(),
                                    ],
                                recombination_rate=recomb_rate, 
                                record_provenance=False, 
                                random_seed=randseed)
    else:
        ts = msprime.sim_ancestry(pops_dict, 
                                sequence_length=seqlen, 
                                demography = m, 
                                recombination_rate=recomb_rate, 
                                record_provenance=False, 
                                random_seed=randseed)

    return ts

# params = one parameter combination.
# mutation rates are entered as a separate list, since they aren't independent coalescent simulations
# n_sims = number of sims to perform per parameter combination
def run_sims(
        params, 
        mut_rates, 
        raster, 
        admix_rast,  
        n_sims, 
        miss_data = 0, 
        transformation = "linear", 
        multi_merge = False, 
        out_prefix = "sims_out.csv"):

    param_df_list = []
    coal_array_list = []
    coal_id_list = []
    demo_summary_df_list = []

    demo_start = time.time()
    # migration rate is Nm / N
    mig_rate = params[1] / params[0]
    demography, demo_summary = setup_demography(
        raster=raster, 
        transformation=transformation, 
        max_k = params[0], 
        admix_rast = admix_rast, 
        merge_time = params[6], 
        ancestral_size_list = params[2], 
        ancestral_merge_size_list=None, 
        ancestral_merge_time_list=params[3], 
        mig_rate = mig_rate, 
        thresh=params[5]
        )
    print("demography setup")
    demo_end = time.time()
    demo_diff = demo_end - demo_start
    seeds = np.random.randint(1, 10000000, size=n_sims)
    # add an id for each parameter combo for easy indexing
    param_id = np.random.randint(1, 1000000)

    for sim in range(n_sims):
        seed = seeds[sim]
        try:
            sim_start = time.time()
            tseq = sim_ts(demography, recomb_rate=params[4], multi_merge=multi_merge, duration=params[6], randseed=seed)
            print("tree sequence simulated")
            sim_end = time.time()
            sim_diff = sim_end - sim_start
            mut_seeds = np.random.randint(1, 1000000, size = len(mut_rates))
            
            tseq = msprime.sim_mutations(tseq, rate = mut_rates[0], random_seed=mut_seeds[0])
            
            coal_list = []
            for j in range(tseq.num_populations):
                # Get the samples corresponding to this population
                samples = tseq.samples(population=j)
                # Simplify the tree sequence to just these samples
                ts_pop = tseq.simplify(samples=samples)
                tree = ts_pop
                # only calc diversity if there is at least two individuals
                if tree.num_samples > 1:
                    coal_list.append(tseq.diversity(samples))
                else: coal_list.append(-1)

            coal_1d = np.array(coal_list)[:-3]

            coal_array = np.reshape(coal_1d, newshape = raster.shape)
            
            print("coalescent times calculated")
            params_df = pd.DataFrame(
                data = {
                "max_k": params[0],
                "Nm": params[1],
                "mig_rate": mig_rate,
                "thresh": params[5],
                "recomb_rate": params[4],
                "mut_rate": mut_rates[0],
                "spatial_merge_time": params[6],
                "demo_time": demo_diff,
                "sim_time": sim_diff,
                "ts_seed": seed,
                "mut_seed": mut_seeds[0],
            }, index = [0]
            )
            

            # add parameter id
            params_df["param_id"] = f"p{param_id}"
            demo_summary["param_id"] = f"p{param_id}"

            param_df_list.append(params_df)
            coal_array_list.append(coal_array)
            coal_id_list.append(param_id)
            demo_summary_df_list.append(demo_summary)
            print("dataframes appended")
                
        except Exception as error:
            print(f"Sim with seed {seed} failed. {error}")
            continue

    # combine parameters into single DataFrame
    if len(param_df_list) > 0:
        param_df_full = pd.concat(param_df_list)
        demo_summary_df_full = pd.concat(demo_summary_df_list)

        print("dataframe lists concatenated")

        # write to files
        out_params = f"{out_prefix}_params.csv"
        out_demo = f"{out_prefix}_demo_summary.csv"
        
        if isfile(out_params):
            param_df_full.to_csv(out_params, mode="a", index=False, header=False)
        else:
            param_df_full.to_csv(out_params, index=False)
        
        print("wrote params to file")

        for i in range(len(coal_array_list)):
            out_coal = f"{out_prefix}_{coal_id_list[i]}_coal.csv"
            np.savetxt(out_coal, coal_array_list[i], delimiter=",")

        if isfile(out_demo):
            demo_summary_df_full.to_csv(out_demo, mode="a", index=False, header=False)
        else:
            demo_summary_df_full.to_csv(out_demo, index=False)
        
        print("write demography summaries to file")
        
        print(f"Finished simulation parameter combo id = {param_id}")
    else:
        print(f"Simulation parameter combo id = {param_id} failed.")

    # remove unnecessary objects to free up memory and not create any weird situations because tseq and demography are mutable
    del tseq
    del demography
    del demo_summary
    
        
#########################################
########## Simulations ##################
#########################################

if __name__=="__main__":
  ##### Read in files and get necessary inputs #####
  raster_path = "projections_cat_inland-removed.tif"
  admix_rast_path = "admix_id_rast_cat.tif"

  with rasterio.open(raster_path) as src:
      r = src.read(masked=True)

  r2 = rasterio.open(raster_path)

  with rasterio.open(admix_rast_path) as src:
      admix_rast = src.read()


  localities_path = "enyalius_locs_genetics.csv"
  localities = pd.read_csv(localities_path)

  localities = localities[localities["species"] == "catenatus"]

  # remove inland individuals
  inland_inds = ["cat_jeq_MTR17108", "cat_jeq_MTR17172", "cat_jeq_MTR17377", "cat_mig_MTR19889", "cat_mig_RPD188", "cat_rem_MTR19915", "cat_rem_RPD128"]
  localities = localities[~localities['id_code'].isin(inland_inds)]

  out_prefix = "sims_2024-04-02_posterior_cat"

  coords = list(zip(localities.longitude, localities.latitude))

  max_thresh = max_thresh_from_coords(r2, coords)

  #admix_dict = admix_to_deme_dict(admix_rast_path=admix_rast_path, pops_dict=pop_dict_simulation)


  # parameter combinations
  # using ranges from posterior estimates for each species. 
  # assuming that the error is normally distributed and approximating the SD as (upper 95 - lower 95) / 4
  n_p = 10000
  
  max_k = np.round(np.random.normal(loc=3083, scale=((4986-1180)/4), size=n_p), 0)

  recomb_rate = [1e-9]
  
  mut_rates = [1e-9] 
  
  # there can't be negative migration, so have to take the absolute value of any negative numbers
  Nm = np.abs(np.random.normal(loc=0.59029424, scale=((2.01201677 + 0.86011314)/4), size=n_p))
  # spatial merge time
  sp_merge_time_list = [22000]

  thresh = [max_thresh]

  ancestral_size_a1 = np.round(np.random.normal(loc=1359495, scale=((2139021-571253)/4), size=n_p), 0)

  # can't have negative population sizes
  ancestral_size_a2 = np.abs(np.round(np.random.normal(loc=503166, scale=((1013851+4341.65625)/4), size=n_p), 0))


  ancestral_size_list = [[ancestral_size_a1[i], ancestral_size_a2[i], ancestral_size_a2[i]] for i in range(len(ancestral_size_a1))]

  #ancestral_merge_size_list = [[ancestral_size_a1[i]] for i in range(len(ancestral_size_a1))]

  ancestral_merge_time_list = [[1000000]]

  ## missing data percentage to correspond with empirical data
  miss_data = 0.4

  # make list of parameter combinations to iterate over
  # lin_sim_params = list(itertools.product(*[
  #     max_k, 
  #     Nm, 
  #     ancestral_size_list, 
  #    #ancestral_merge_size_list, 
  #     ancestral_merge_time_list, 
  #     recomb_rate, 
  #     [0],
  #     sp_merge_time_list]))
  num_samples = 500
  hinge_sim_params = list(zip(np.random.choice(max_k, num_samples),
                              np.random.choice(Nm, num_samples),
                              random.choices(ancestral_size_list, k=num_samples),
                              random.choices(ancestral_merge_time_list, k=num_samples),
                              np.random.choice(recomb_rate, num_samples),
                              np.random.choice(thresh, num_samples),
                              np.random.choice(sp_merge_time_list, num_samples)))
  

  # number of coalescent simulations per parameter combination
  n_sims_hinge = 1
  # since hinge has an extra parameter that varies (threshold), multiply the number of hinge simulations per parameter combo by the number of threshold parameters explored
  #n_sims_lin = n_thresh * n_sims_hinge
  n_sims_lin = 1

  # iterate over simulations

  
  run_sim_part = partial(run_sims,
                          mut_rates=mut_rates,
                          raster=r,
                          admix_rast=admix_rast,
                          # merge_time = 22000,
                          n_sims=n_sims_hinge,
                          miss_data=miss_data,
                          transformation="hinge",
                          multi_merge=False,
                          out_prefix=out_prefix)

  pool = Pool(20)
  
  pool.map(run_sim_part, hinge_sim_params)

  pool.close()
  pool.join()
  
