###### Exploratory simulations to see impacts of parameter values on summary statistics #######

# import packages
from os.path import isfile
import numpy as np
from numpy import ma
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
import numba
import tracemalloc
from mpi4py import MPI

# initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#### Define functions #####

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
def stepping_stone2d(initial_size, rate, scale=True):
    assert len(initial_size.shape) <= 3

    n, m = initial_size.shape
    N = n * m
    model = msprime.Demography.isolated_model(initial_size.reshape(N))

    # set population names
    for j in range(n):
        for k in range(m):
            index = j * m + k
            model.populations[index].name = f"pop_{j}_{k}"
    
    # setup migration rate matrices
    if np.array(rate).ndim == 0:
        if scale:
            model.migration_matrix = migration_matrix(initial_size, rate, scale=True)
        else: 
            model.migration_matrix = migration_matrix(initial_size, rate, scale=False)
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

## add_nonspatial_phase ##
def add_nonspatial_phase(model, ancestral_size, merge_time):
    # add an ancestral population
    model.add_population(name = "ANC", initial_size=ancestral_size)

    # get names of populations for initiating the collecting phase
    pop_names = [[pop.name] for pop in model.populations if pop.name != "ANC"]
    
    # add the time when the spatial simulation collapses into the collecting phase
    [model.add_population_split(time = merge_time, derived = name, ancestral = "ANC") for name in pop_names]

    return model


## coords_to_sample_dict ##
### In this case, coordinates are a list of Point geometries (in the coords_to_sample_dict_empirical they are a list of tuples. I NEED TO COMBINE THE TWO FUNCTIONS SO IT ISN'T CONFUSING)
def coords_to_sample_dict(raster, coordinates, num_inds):

    # mask the raster, so only cells that the coordinates overlap are masked
    out_image = mask(raster, coordinates, nodata="nan", filled = False)

    # unravel the masked raster into a single dimension, only keeping the first raster in case there are multiple rasters in the dataset
    oi_1d = out_image[0][0].ravel()

    # get the indices of the localities
    cell_id = ma.nonzero(oi_1d)[0]
    

    if isinstance(num_inds, int):

        # function to split list into even-sized chunks from a list
        # stolen from https://stackoverflow.com/a/312464/4406870
        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), n):
                yield lst[i:i + n]

        # short pop dictionary for simulation sampling
        pop_dict = {key: num_inds for key in cell_id}
        
        # make a longer pop dictionary for calculating allele counts per subpopulation
        total_inds_ind = list(range(num_inds * len(cell_id)))

        # split into even-sized chunks
        inds_split = chunks(total_inds_ind, num_inds)

        # make into a dictionary
        pop_dict_long = dict(zip(cell_id, inds_split))

        pop_dict_long


    
    return pop_dict, pop_dict_long

# coordinates should be a list of longitude, latitude tuples corresponding with each individual in your genetic dataset
# individual_ids should be the unique identification number of each individual at each coordinate in `coordinates`. The individual_id should match its cid in the vcf (i.e. the "samples" slot in you vcf should link up with the individual_ids)
# outputs a "pop_dict_simulation" object, which outputs the number of individuals to sample per population for the simulations to align with the empirical sampling
def coords_to_sample_dict_empirical(raster, coordinates):
    # get the cell that each individual belongs to  
    ## I have to iterate through each locality to get the cell ID for each individual locality. Masking everything returns the cell IDs out of order.
    cell_list = []
    for xy in coordinates:
    
        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]
        
        # mask the raster with the points
        out_image = mask(raster, pt2, nodata="nan", filled = False)

        # turn into 1D array
        oi_1d = out_image[0][0].ravel()

        # # get the indices of the locality and append
        cell_id = ma.nonzero(oi_1d)[0]    
        cell_list.append(cell_id)
    
    cell_id_array = np.concatenate(cell_list, axis = 0)

    # get the number of individuals to sample from the simulation
    pop_dict_simulation = {}
    for cid in np.unique(cell_id_array):
        num_inds = np.sum(cell_id_array == cid) 
        pop_dict_simulation[cid] = num_inds

    pop_dict_sim_long = {}
    unique_ids = np.unique(cell_id_array)
    for i in range(len(unique_ids)):
        cid = unique_ids[i]
        if i > 0:
            cid_prev = unique_ids[i - 1]
            last_ind_prev = pop_dict_sim_long[cid_prev][len(pop_dict_sim_long[cid_prev]) - 1]
            pop_dict_sim_long[cid] = np.array(range(last_ind_prev + 1, last_ind_prev + pop_dict_simulation[cid] + 1))
        else: 
            pop_dict_sim_long[cid] = np.array(range(pop_dict_simulation[cid]))
    
    return pop_dict_simulation, pop_dict_sim_long

## max_thresh_from_coords ##
def max_thresh_from_coords(raster, coordinates):

    # mask the raster, so only cells that the coordinates overlap are masked
    out_image = mask(raster, coordinates, nodata="nan", filled = False)

    # find the minimum value of the masked raster
    max_thresh = np.min(out_image[0])

    return max_thresh


## counts_from_ts ##
def counts_from_ts(ts, pops_dict, missing_data_perc = 0, r2_thresh = 0.1, seed = 1): 
    np.random.seed(seed)
    # create data mask that masks the tree sequence genotype matrix to have similar missing data patterns to the empirical data

    # get total the number of elements in the genotype matrix
    total_array_size = ts.genotype_matrix().shape[0] * ts.genotype_matrix().shape[1]

    # fill an array with zeros that is the total size of the original matrix
    a = np.zeros(total_array_size, dtype=int)

    # set the first x elements to be the size of non-missing data
    non_miss = int(np.ceil(total_array_size * (1 - missing_data_perc)))

    # 1 means "True"
    a[:non_miss] = 1

    # randomly shuffle True and False values
    np.random.shuffle(a)

    # transform to boolean
    a = a.astype(bool)

    # reshape to genotype matrix shape
    miss_mask = np.reshape(a, ts.genotype_matrix().shape)

    # mask the genotype matrix
    geno_mat = np.ma.masked_array(ts.genotype_matrix(), mask=miss_mask, fill_value=-1)

    # convert to a GenotypeArray
    h = allel.HaplotypeArray(geno_mat)

    gt = h.to_genotypes(ploidy = 2)

    # get the allele counts matrix for all individuals
    gc = gt.count_alleles(max_allele = 1)

    # only retain segregating alleles
    is_seg = gc.is_segregating()
    gt_seg = gt.compress(is_seg, axis=0)
    gc_seg = gc.compress(is_seg)

    # alt allele counts matrix
    # needed for finding unlinked SNPs
    gn = gt_seg.to_n_alt(fill=-1)

    # filter for unlinked loci
    # locate unlinked SNPs
    loc_unlinked = allel.locate_unlinked(gn, threshold=r2_thresh)

    # select unlinked SNPs
    gc_unlinked = gc_seg[loc_unlinked]
    gt_unlinked = gt_seg[loc_unlinked]

    # remove singletons
    not_singleton = ~gc_unlinked.is_singleton(allele = 1)
    gc_unlinked = gc_unlinked[not_singleton]

    gc_unlinked_pops = None

    if pops_dict:
        gt_unlinked = gt_unlinked[not_singleton]
        gc_unlinked_pops = gt_unlinked.count_alleles_subpops(max_allele = 1, subpops = pops_dict)

    return gc_unlinked, gc_unlinked_pops

## hill_number ##
def hill_number(freqs, order):
    if order == 0:
        return len(np.nonzero(freqs)[0])
    if order == 1:
        h1 = np.exp(stats.entropy(freqs))
        return h1
    tot = float(np.sum(freqs))
    proportions = np.array(freqs[freqs > 0])/tot
    prop_order = proportions**order
    h2 = np.sum(prop_order)**(1/(1-order))
    return h2

## sampled_cells_to_coords ##
def sampled_cells_to_coords(raster, coordinates):
    
    # mask the raster, so only cells that the coordinates overlap are masked
    out_image = mask(raster, coordinates, nodata="nan", filled = False)

    # get the indices of the localities
    cell_ind = [tuple(_) for _ in np.transpose(ma.nonzero(out_image[0][0]))]


    rows, cols = zip(*cell_ind)  # unzip the cell indices to separate row and col lists
    xs, ys = rasterio.transform.xy(raster.transform, rows, cols)
    lons = np.array(xs)
    lats = np.array(ys)

    return list(zip(lons, lats))


## calc_sumstats ##

def calc_sumstats(counts_array, counts_array_pops, pops_coordinates):

    sfs = allel.sfs_folded(counts_array)

    # calculate the first 3 Hill numbers of the site frequency spectrum, scaling by the sample size
    sfs_h1 = hill_number(sfs, 1) / len(sfs)
    sfs_h2 = hill_number(sfs, 2) / len(sfs)
    sfs_h3 = hill_number(sfs, 3) / len(sfs)


    # number of variants
    num_var = counts_array.shape[0]

    # average pi across sites
    pi = np.nanmean(allel.mean_pairwise_difference(counts_array, fill = np.nan))


    # tajima's D 
    taj_d = np.nanmean(allel.moving_tajima_d(counts_array, size=100, step=10))

    stat_df = pd.DataFrame(
        data = {
            "sfs_h1": sfs_h1,
            "sfs_h2": sfs_h2,
            "sfs_h3": sfs_h3,
            "num_var": num_var,
            "pi": pi,
            "taj_d": taj_d,
        }, index = [0]
    )

    
    if counts_array_pops:

        # get pi for all pops
        pi_pops = []
        for key in counts_array_pops:
            pi_pops.append(np.nanmean(allel.mean_pairwise_difference(counts_array_pops[key])))
        
        pi_pops_dict = dict(zip(counts_array_pops.keys(), pi_pops))

        # add pi to dataframe
        for key in pi_pops_dict:
            colname = f"pi_pop_{key}"
            stat_df[colname] = pi_pops_dict[key]

        #################################
        #### ISOLATION BY DISTANCE ######
        #################################

        ##### slope and r2 of gen ~ geo regression ##### (similar to mantel corr)
        # get pairwise dxy for all pops
        dxy = []
        for ac_ind in list(itertools.combinations(list(counts_array_pops.keys()), 2)):
            ac1 = counts_array_pops[ac_ind[0]]
            ac2 = counts_array_pops[ac_ind[1]]
            
            d = np.nanmean(allel.mean_pairwise_difference_between(ac1, ac2))
            dxy.append(d)

        # scale dxy according to Roussett 1997 (they used FST, but logic still follows)
        dxy = np.asarray(dxy)
        
        dxy_scaled = dxy / (1 - dxy)
        
        # get pairwise geographic distance for all pops
        
        # convert to radians
        long_rad = [radians(x[0]) for x in pops_coordinates]
        lat_rad = [radians(y[1]) for y in pops_coordinates]

        geometry = list(zip(long_rad, lat_rad))
        geometry = pd.unique(geometry)
        geometry = [list(_) for _ in geometry]

        dist = haversine_distances(geometry) * 6371000/1000  # multiply by Earth radius to get kilometers

        # get the upper diagonal of the distance matrix
        dist = dist[np.triu_indices(dist.shape[0], k=1)]

        
        # take the log10 of geographic distance to scale linearly for IBD regression
        logdist = np.log10(dist)

        dist_df = pd.DataFrame(
            data = {
                "geo_dist": logdist,
                "dxy": dxy_scaled
            }
        )

        dist_df = dist_df.dropna()

        # linear model, extracting the R2 and coefficient (slope)
        reg = linear_model.LinearRegression()

        # reshape X so it is 2D
        geo_dist = np.array(dist_df["geo_dist"])
        geo_dist = geo_dist.reshape(-1, 1)
        
        reg.fit(geo_dist, dist_df["dxy"])

        r2 = reg.score(geo_dist, dist_df["dxy"])
        b = reg.coef_
        
        stat_df["ibd_r2"] = r2
        stat_df["ibd_slope"] = b

        ####### Moran's I ########

        ## geographic weights from Voronoi polygons
        coords_rad = np.array([long_rad, lat_rad]).transpose()

        weights = Voronoi(coords_rad, use_index=False)
        # array of per-pop pi, calculated earlier
        pi_array = np.array(pi_pops)
        moran_i = Moran(pi_array, weights)
        
        mi_i = moran_i.I

        stat_df["morans_i"] = mi_i
        

    return stat_df


# setup demography
# returns a tuple of a Demography object and a summary of population sizes
def setup_demography(raster, transformation, max_k, ancestral_n, ancestral_n_multiplier, mig_rate = 0.001, thresh=0, tstep = 1000):

    k = raster_to_k(raster=raster, transformation=transformation, threshold = thresh, max_local_k=max_k)

    m = stepping_stone2d(k[0], rate = mig_rate)

    m = add_landscape_change(model = m, k_stack = k, rate=mig_rate, timestep = tstep)

    # if there is no argument supplied to the ancestral_n argument, make the ancestral_n equal to the total number of individual across the landscape at the most ancient time period
    if ancestral_n is None and ancestral_n_multiplier is not None:
        ancestral_n = np.nansum(k[k.shape[0] - 1]) * ancestral_n_multiplier
        m = add_nonspatial_phase(m, ancestral_size=ancestral_n, merge_time=k.shape[0] * tstep + tstep)
    elif ancestral_n is None and ancestral_n_multiplier is None: 
        ancestral_n = np.nansum(k[k.shape[0] - 1])
        m = add_nonspatial_phase(m, ancestral_size=ancestral_n, merge_time=k.shape[0] * tstep + tstep)
    else:
        m = add_nonspatial_phase(m, ancestral_size=ancestral_n, merge_time=k.shape[0] * tstep + tstep)

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
        avg_nonzero_deme_size = np.nanmean(tstep_nz)
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
    total_inds_avg = np.mean(total_inds_list)
    total_inds_sd = np.std(total_inds_list)

    ## average and standard deviation of average number of average number of individuals per deme
    avg_deme_size_avg = np.mean(deme_size_list)
    avg_deme_size_sd = np.std(deme_size_list)

    ## number of occupied demes in the current landscape
    num_occ_curr = num_occupied_list[0]
    ## number of occupied demes in the landscape with the largest number of occupied demes and its index
    num_occ_largest = np.max(num_occupied_list)
    num_occ_largest_index = np.argmax(num_occupied_list)
    ## average number of occupied demes across timesteps
    num_occ_average = np.mean(num_occupied_list)
    ## standard deviation of the number of occupied demes
    num_occ_sd = np.std(num_occupied_list)

    # pandas DataFrame of the demography summaries
    demo_df = pd.DataFrame(
            data = {
                "ancestral_n": ancestral_n,
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
    

    return m, demo_df


# m is the demography
def sim_ts(m, recomb_rate, mut_rate, pops_dict="none", randseed=1):
# if there is no user-supplied pops_dict, sample all populations with at least two individuals
    if pops_dict == "none":
        # get populations with at least two individuals to sample for coalescent simulations
        nzp = []
        
        for i in range(len(m.populations)):
            if m.populations[i].initial_size >= 2.0 and m.populations[i].name != "ANC":
                nzp.append(m.populations[i].name)

        pops_dict = {key: 2 for key in nzp}

    ts = msprime.sim_ancestry(pops_dict, sequence_length=1e6, demography = m, recombination_rate=recomb_rate, record_provenance=False, random_seed=randseed)

    ts = msprime.sim_mutations(ts, rate = mut_rate, random_seed=randseed)

    return ts

# params = one parameter combination.
# n_sims = number of sims to perform per parameter combination
def run_sims(params, raster, pops_dict, pops_dict_long, pops_coordinates, n_sims, miss_data = 0, transformation = "linear", out_file = "sims_out.csv"):

    df_list = []
    demo_start = time.time()
    # migration rate is Nm / N
    mig_rate = params[1] / params[0]
    demography, demo_summary = setup_demography(raster=raster, transformation=transformation, max_k = params[0], ancestral_n = None, ancestral_n_multiplier = params[2], mig_rate = mig_rate, thresh=params[5])
    demo_end = time.time()
    demo_diff = demo_end - demo_start
    seeds = np.random.randint(1, 10000000, size=n_sims)

    for sim in range(n_sims):
        seed = seeds[sim]
        try:
            sim_start = time.time()
            tseq = sim_ts(demography, recomb_rate=params[3], mut_rate=params[4], pops_dict=pops_dict, randseed=seed)
            sim_end = time.time()
            sim_diff = sim_end - sim_start
            c, c_long = counts_from_ts(ts=tseq, pops_dict=pops_dict_long, missing_data_perc=miss_data, seed=seed)
            sumstats = calc_sumstats(c, c_long, pops_coordinates = pops_coordinates)
            params_df = pd.DataFrame(
                data = {
                "max_k": params[0],
                "Nm": params[1],
                "ancestral_n_multiplier": params[2],
                "mig_rate": mig_rate,
                "thresh": params[5],
                "recomb_rate": params[3],
                "mut_rate": params[4],
                "demo_time": demo_diff,
                "sim_time": sim_diff,
                "seed": seed
            }, index = [0]
            )
            
            full_df = pd.concat([params_df, sumstats], axis=1)
            df_list.append(full_df)
            del tseq
        except Exception as error:
            print(f"Sim with seed {seed} failed. {error}")
    
    # combine parameters into single DataFrame
    df_full = pd.concat(df_list)

    # add an id for each parameter combo for easy indexing
    param_id = np.random.randint(1, 1000000)
    df_full["param_id"] = f"p{param_id}"

    # add demographic summaries to df
    df_full.reset_index(inplace=True, drop=True) # need to reset indices so I can concatenate the two DataFrames
    df_full = pd.concat([df_full, demo_summary], axis=1)

    # remove demography object
    del demography
    del demo_summary
    
    # write to file
    if isfile(out_file):
        df_full.to_csv(out_file, mode="a", index=False, header=False)
    else:
        df_full.to_csv(out_file, index=False)
    
    print(f"Finished simulation parameter combo id = {param_id}")

#########################################
########## Simulations ##################
#########################################

def main():
    ##### Read in files and get necessary inputs #####
    raster_path = "projections_ihe.tif"

    with rasterio.open(raster_path) as src:
        r = src.read(masked=True)

    r2 = rasterio.open(raster_path)

    localities_path = "enyalius_locs_genetics.csv"
    localities = pd.read_csv(localities_path)

    localities = localities[localities["species"] == "iheringii"]

    out_file_linear = "sims_2023-11-22_linear_second.csv"
    out_file_hinge = "sims_2023-11-18_hinge.csv"

    geometry = [Point(xy) for xy in zip(localities.longitude, localities.latitude)]

    pop_coords = sampled_cells_to_coords(r2, geometry)

    pop_dict_simulation, pop_dict_sim_long = coords_to_sample_dict_empirical(r2, zip(localities.longitude, localities.latitude))

    max_thresh = max_thresh_from_coords(r2, geometry)


    # parameter combinations
    n_p = 10
    n_thresh = 3
    max_k = np.linspace(100, 1000, n_p,  dtype = int)
    recomb_rate = [1e-9]
    mut_rate = [1e-8]
    Nm = np.geomspace(0.1, 10, n_p)
    thresh = np.linspace(max_thresh / 4, max_thresh, n_thresh)
    ## I'm setting the ancestral N value as fractions of the total number of individuals across the landscape
    ## This accounts for population structure/varying migration rates impacting the ancestral Ne
    ancestral_n_multiplier = [0.01, 0.1, 1]

    ## missing data percentage to correspond with empirical data
    miss_data = 0.4

    # make list of parameter combinations to iterate over
    lin_sim_params = list(itertools.product(*[max_k, Nm, ancestral_n_multiplier, recomb_rate, mut_rate, [0]]))
    hinge_sim_params = list(itertools.product(*[max_k, Nm, ancestral_n_multiplier, recomb_rate, mut_rate, thresh]))

    # number of simulations per parameter combination
    n_sims_hinge = 5
    # since hinge has an extra parameter that varies (threshold), multiply the number of hinge simulations per parameter combo by the number of threshold parameters explored
    n_sims_lin = n_thresh * n_sims_hinge

    # iterate over simulations

    ## I'm separating them into batches because for some reason the Huxley cluster doesn't want to split my jobs across nodes
    mid_ind = len(lin_sim_params) // 2
    lin_sim_params_first = [lin_sim_params[i] for i in range(mid_ind + 1)]
    lin_sim_params_second = [lin_sim_params[i] for i in range(mid_ind, len(lin_sim_params))]

    ### linear sims
    num_param_combos = len(lin_sim_params_second)
    sims_per_process = num_param_combos // size
    if num_param_combos % size > 0:
        sims_per_process += 1
    
    first_sim = rank * sims_per_process
    last_sim = first_sim + sims_per_process

    for sim in range(first_sim, last_sim):
        run_sims(params = lin_sim_params_second[sim], raster = r, pops_dict = pop_dict_simulation, pops_dict_long = pop_dict_sim_long, pops_coordinates = pop_coords, n_sims = n_sims_lin, miss_data = miss_data, out_file = out_file_linear)

if __name__=="__main__":

    tracemalloc.start()
    start_sims_time = time.time()
    main()
    end_sims_time = time.time()
    total_time_hours = (end_sims_time - start_sims_time) / 3600
    current_mem, max_mem = tracemalloc.get_traced_memory()
    print(f"{current_mem=}, {max_mem=}")
    print(f"Total time = {total_time_hours} hours")