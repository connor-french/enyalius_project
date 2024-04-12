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
from scipy.interpolate import NearestNDInterpolator
import numba
import tracemalloc
from scipy.stats import entropy
import warnings
from mpi4py import MPI

# ignoring the performance warnings from how I build the pandas dataframes
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

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
        avg_sdm = np.nanmean(t[t>=threshold])
        t[t >= threshold] = np.ceil(avg_sdm * max_local_k)
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
# raster is a raster opened in rasterio in "r" mode and is the same raster used to parameterize the demographic model
# coordinates is a shapely Point geometry that can be obtained using the shapely.geometry module or as a geometry column in a GeoPandas dataframe. Note- make sure the raster and coordinates are in the same projection.
# num_inds is either a single number or a vector of numbers that correspond with the number of individuals (not genomes) you want to sample per population

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
def coords_to_sample_dict_empirical(raster, coordinates, individual_ids=None, vcf_path=None):
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
    
    if individual_ids is None:
        try:
            if vcf_path is not None:
                raise ValueError
        except ValueError:
            print("When a VCF path is provided, individual IDs corresponding to those in the VCF are expected. Please provide those IDs.")
    elif vcf_path is None:
        try:
            if individual_ids is not None:
                raise ValueError
        except ValueError:
            print("When individual IDs are provided, a VCF file is expected. Please provide a path to your VCF file.")
    else:
        pop_dict_empirical = None
    
    if individual_ids is not None and vcf_path is not None:
        
        # read in the vcf
        callset = allel.read_vcf(vcf_path)

        # get the indices of each individual in the vcf
        ind_indices = []
        for cid in individual_ids:
            index = list(callset["samples"]).index(cid) 
            ind_indices.append(index)

        ind_indices_array = np.array(ind_indices)
        

        # make each cell id a key and the individuals that belong to that cell id the values
        pop_dict_empirical = {}
        for cid in np.unique(cell_id_array):
            id_inds = np.where(cell_id_array ==cid) 
            id_ind_indices = ind_indices_array[id_inds]
            pop_dict_empirical[cid] = id_ind_indices
    else:
        pop_dict_empirical = None
    
    return pop_dict_simulation, pop_dict_sim_long, pop_dict_empirical

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


## counts_from_ts ##
# ts = the msprime tree sequence from the simulation
# pops_dict_deme = a dictionary linking sampled individuals to their deme
# pops_dict_admix = a dictionary linking deme IDs to their admixture population assignment
# missing_data_perc = the percentage of missing data observed in empirical data
# r2_thresh = the linkage disequilibrium threshold to use for retaining unlinked SNPs

# returns 
def counts_from_ts(ts, pops_dict_deme=None, pops_dict_admix=None, missing_data_perc = 0, r2_thresh = 0.1, seed = 1): 
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
    gt_unlinked = gt_unlinked[not_singleton]

    if pops_dict_deme is not None:
        gc_unlinked_demes = gt_unlinked.count_alleles_subpops(max_allele = 1, subpops = pops_dict_deme)
    else:
        gc_unlinked_demes = None
    

    # create allele counts matrices for the admixture populations
    if pops_dict_admix is not None:
        admix_dict_inds = {}
        for key in pops_dict_admix:
            demes = pops_dict_admix[key]
            admix_dict_inds[key] = []
            for d in demes:
                # convert numpy array to list
                inds = pops_dict_deme[d].tolist()
                admix_dict_inds[key].extend(inds)

        gc_unlinked_admix = gt_unlinked.count_alleles_subpops(max_allele = 1, subpops = admix_dict_inds)
        
    else:
        gc_unlinked_admix = None    


    return gc_unlinked, gc_unlinked_demes, gc_unlinked_admix

## sampled_cells_to_coords ##
# raster is a raster read into rasterio in "r" mode
# coordinates are a list of coordinate tuples
# returns a dictionary, where the keys are the population IDs and the values are a list of the long-lat coordinate pairs
def sampled_cells_to_coords(raster, coordinates):
    
    # get the population IDs
    ## I have to iterate through each locality to get the cell ID for each individual locality. Masking everything returns the cell IDs out of order.
    ## I'm creating a dictionary 
    cell_dict = {}
    for xy in coordinates:

        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]
        
        # mask the raster with the points
        out_image = mask(raster, pt2, nodata="nan", filled = False)

        ###### GET CELL INDICES ######
        # turn into 1D array
        oi_1d = out_image[0][0].ravel()

        # # get the indices of the locality and append
        cell_id = ma.nonzero(oi_1d)[0][0] 

        ####### GET CELL COORDINATES #############
        # get the indicex of the locality
        cell_ind = [tuple(_) for _ in np.transpose(ma.nonzero(out_image[0][0]))]

        rows, cols = zip(*cell_ind)  # unzip the cell indices to separate row and col lists
        coords_tuple = rasterio.transform.xy(raster.transform, rows, cols)
        
        # I have to convert the x and y coords to a list instead of a tuple, because the transform.xy function returns a tuple where each long and lat is a single-value list, which makes indexing later convoluted.
        x = coords_tuple[0][0]
        y = coords_tuple[1][0]

        coords = [x, y]

        cell_dict[cell_id] = coords
        
    return cell_dict

## calc_sumstats ##


# counts_array is an allele counts array of the entire simulated population.
# demes_coordinates is a a dictionary, where the deme ID is the key and a list of the lat-long coordinate associated is the value.
# between_admix_pop_sumstats is a boolean that indicates whether you want to calculate between-population summary statistics like FST and DXY for the admixture populations. The default is False, because we're assuming the user is interested in the more recent, spatially-explicit period, rather than more ancient genetic structure.

def calc_sumstats(counts_array, demes_coordinates, admix_dict=None, counts_array_demes=None, counts_array_admix=None, between_admix_pop_sumstats = False):
    
    # calculate scaled Hill numbers of the SFS, using Isaac Overcast's approach
    def hill_number(freqs, order):
        if order == 0:
            return len(np.nonzero(freqs)[0])
        if order == 1:
            h1 = np.exp(entropy(freqs))
            return h1
        tot = float(np.sum(freqs))
        proportions = np.array(freqs[freqs > 0])/tot
        prop_order = proportions**order
        h2 = np.sum(prop_order)**(1/(1-order))
        return h2
    
    # get pairwise dxy for pops
    # deme IDs is a list of population IDs we want to calculate DXY for.
    def calc_dxy(deme_ids):
        dxy = []
        ## also get the names of each pop pair to identify the raw dxy
        dxy_name = []

        # I'm using the demes_coordinates dictionary to iterate through the population names, so I can be sure that the DXY distances match up with the geographic distances
        for ac_ind in list(itertools.combinations(deme_ids, 2)):
            ac1 = counts_array_demes[ac_ind[0]]
            ac2 = counts_array_demes[ac_ind[1]]
            
            d = np.nanmean(allel.mean_pairwise_difference_between(ac1, ac2))
            d_name = f"dxy_{ac_ind[0]}_{ac_ind[1]}"

            dxy.append(d)
            dxy_name.append(d_name)
        return dxy, dxy_name

    # isolation-by-distance function
    ##### slope and r2 of gen ~ geo regression ##### (similar to mantel corr)
    def calc_ibd(dxy, coords):
        
        # scale dxy according to Roussett 1997 (they used FST, but logic still follows)
        dxy = np.asarray(dxy)
        
        dxy_scaled = dxy / (1 - dxy)
        
        # get pairwise geographic distance for all pops
        
        # convert to radians
        long_rad = [radians(x[0]) for x in coords]
        lat_rad = [radians(y[1]) for y in coords]

        geometry = list(zip(long_rad, lat_rad))

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

        return b, r2


    
    # only calculate species-wide stats if there aren't admixture populations
    if counts_array_admix is None:
        sfs = allel.sfs_folded(counts_array)

        # calculate the first 5 Hill numbers of the site frequency spectrum, scaling by the sample size
        sfs_h1 = hill_number(sfs, 1) / len(sfs)
        sfs_h2 = hill_number(sfs, 2) / len(sfs)
        sfs_h3 = hill_number(sfs, 3) / len(sfs)
            
        # average pi across sites
        pi_raw = allel.mean_pairwise_difference(counts_array, fill = np.nan)
        pi = np.nanmean(pi_raw)
        # standard deviation of pi across sites
        pi_sd = np.nanstd(pi_raw)


        # tajima's D 
        taj_d_raw = allel.moving_tajima_d(counts_array, size=100, step=10)
        taj_d = np.nanmean(taj_d_raw)
        taj_d_sd = np.nanstd(taj_d_raw)

        stat_df = pd.DataFrame(
            data = {
                "sfs_h1": sfs_h1,
                "sfs_h2": sfs_h2,
                "sfs_h3": sfs_h3,
                "pi": pi,
                "pi_sd": pi_sd,
                "taj_d": taj_d,
                "taj_d_sd": taj_d_sd
            }, index = [0]
        )

        # calculate dxy
        dids = list(demes_coordinates.keys())

        dxy, dxy_name = calc_dxy(dids)
        
        for name in range(len(dxy_name)):
            stat_df[dxy_name[name]] = dxy[name] 


    # otherwise, calculate stats per admixture population
    # also calculate Hudson's FST and Dxy between the two populations
    else:
        stat_df = pd.DataFrame(
            data = {},
            index = [0]
        )

        # calc sfs Hill numbers, pi, and tajima's D
        for key in counts_array_admix:
            sfs = allel.sfs_folded(counts_array_admix[key])
            stat_df[f"sfs_h1_a{key}"] = hill_number(sfs, 1) / len(sfs)
            stat_df[f"sfs_h2_a{key}"] = hill_number(sfs, 2) / len(sfs)
            stat_df[f"sfs_h3_a{key}"] = hill_number(sfs, 3) / len(sfs)
            pi_raw = allel.mean_pairwise_difference(counts_array_admix[key])
            stat_df[f"pi_a{key}"] = np.nanmean(pi_raw)
            stat_df[f"pi_sd_a{key}"] = np.nanstd(pi_raw)
            taj_d_raw = allel.moving_tajima_d(counts_array_admix[key], size=100, step=10)
            stat_df[f"taj_d_a{key}"] = np.nanmean(taj_d_raw)
            stat_df[f"taj_d_sd_a{key}"] = np.nanstd(taj_d_raw)
        
        # Dxy and Hudson's FST for admixture populations and Dxy for all pops.
        if between_admix_pop_sumstats:
            for ac_ind in list(itertools.combinations(list(counts_array_admix.keys()), 2)):
                ac1 = counts_array_admix[ac_ind[0]]
                ac2 = counts_array_admix[ac_ind[1]]
                d_raw = allel.mean_pairwise_difference_between(ac1, ac2)
                d = np.nanmean(d_raw)
                d_sd = np.nanstd(d_raw)
                stat_df[f"dxy_a{ac_ind[0]}_{ac_ind[1]}"] = d
                stat_df[f"dxy_sd_a{ac_ind[0]}_{ac_ind[1]}"] = d_sd
                fst, fst_se, _, _ = allel.average_hudson_fst(ac1, ac2, 100)
                stat_df[f"fst_a{ac_ind[0]}_{ac_ind[1]}"] = fst
                stat_df[f"fst_se_a{ac_ind[0]}_{ac_ind[1]}"] = fst_se   
            
            # calculate dxy
            dids = list(demes_coordinates.keys())

            dxy, dxy_name = calc_dxy(dids)
            
            for name in range(len(dxy_name)):
                stat_df[dxy_name[name]] = dxy[name]       
    

    if counts_array_demes is not None:
        # get pi for all pops
        pi_pops = []
        for key in counts_array_demes:
            pi_pops.append(np.nanmean(allel.mean_pairwise_difference(counts_array_demes[key])))
        
        pi_pops_dict = dict(zip(counts_array_demes.keys(), pi_pops))

        # add pi to dataframe
        for key in pi_pops_dict:
            colname = f"pi_pop_{key}"
            stat_df[colname] = pi_pops_dict[key]

        
        if admix_dict is None:

            # convert the values of the demes_coordinates dictionary into a list
            coords = list(demes_coordinates.values())
            b, r2 = calc_ibd(dxy, coords)
            stat_df["ibd_slope"] = b
            stat_df["ibd_r2"] = r2
            
        else:
            
            for key in admix_dict:
                demes = admix_dict[key]

                # get the coordinates associated with the demes
                coords_d = []
                for d in demes:
                    c = demes_coordinates[d]
                    coords_d.append(c)
                
                ## calculate dxy. This is super fast, so it's easier to do this than try weird subsetting to subset the full IBD matrix for the DXY values we want per admix pop
                dxy, dxy_name = calc_dxy(demes)

                for name in range(len(dxy_name)):
                    stat_df[dxy_name[name]] = dxy[name] 

                ## calculate ibd
                b, r2 = calc_ibd(dxy, coords_d)
                stat_df[f"ibd_slope_a{key}"] = b
                stat_df[f"ibd_r2_a{key}"] = r2
    

    ####### Moran's I ########
    # convert the values of the demes_coordinates dictionary into a list
    coords = list(demes_coordinates.values())
    # convert to radians
    long_rad = [radians(x[0]) for x in coords]
    lat_rad = [radians(y[1]) for y in coords]
    ## geographic weights from Voronoi polygons
    coords_rad = np.array([long_rad, lat_rad]).transpose()

    weights = Voronoi(coords_rad, use_index=False)
    
    # array of per-pop pi, calculated earlier
    pi_array = np.array(pi_pops)
    moran_i = Moran(pi_array, weights)
    
    mi_i = moran_i.I

    stat_df["morans_i"] = mi_i

    return np.round(stat_df, 6)

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
    
    demo_df = np.round(demo_df, 6)
    
    for i in range(len(ancestral_size_list)):
        demo_df[f"admix_n_a{i + 1}"] = ancestral_size_list[i]
    
    for i in range(len(ancestral_merge_size_list)):
        demo_df[f"admix_anc_n_a{i + 1}"] = ancestral_merge_size_list[i]
    
    for i in range(len(ancestral_merge_time_list)):
        demo_df[f"anc_merge_time_a{i + 1}"] = ancestral_merge_time_list[i]

    return m, demo_df


# m is the demography
def sim_ts(m, recomb_rate, seqlen = 1e6, pops_dict=None, multi_merge = False, duration = None, randseed=1):
# if there is no user-supplied pops_dict, sample all populations with at least two individuals
    if pops_dict is None:
        # get populations with at least two individuals to sample for coalescent simulations
        nzp = []
        
        for i in range(len(m.populations)):
            if m.populations[i].initial_size >= 2.0 and m.populations[i].name != "ANC":
                nzp.append(m.populations[i].name)

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
def run_sims(params, mut_rates, raster, admix_rast,  pops_dict_sim, pops_dict_deme, pops_dict_admix, demes_coordinates, n_sims, miss_data = 0, transformation = "linear", multi_merge = False, out_prefix = "sims_out.csv"):

    param_df_list = []
    sumstats_df_list = []
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
            tseq = sim_ts(demography, recomb_rate=params[4], pops_dict=pops_dict_sim, multi_merge=multi_merge, duration=params[6], randseed=seed)
            print("tree sequence simulated")
            sim_end = time.time()
            sim_diff = sim_end - sim_start
            mut_seeds = np.random.randint(1, 1000000, size = len(mut_rates))
            sim_counter = 0
            for m in range(len(mut_rates)):
                tseq = msprime.sim_mutations(tseq, rate = mut_rates[m], random_seed=mut_seeds[m])
                print("mutations simulated")
                c, c_demes, c_admix = counts_from_ts(ts=tseq, pops_dict_deme=pops_dict_deme, pops_dict_admix=pops_dict_admix, missing_data_perc=miss_data, seed=mut_seeds[m])
                print("counts arrays created")
                sumstats = calc_sumstats(counts_array=c, demes_coordinates=demes_coordinates, admix_dict=pops_dict_admix, counts_array_demes=c_demes, counts_array_admix=c_admix)
                print("sumstats calculated")
                params_df = pd.DataFrame(
                    data = {
                    "max_k": params[0],
                    "Nm": params[1],
                    "mig_rate": mig_rate,
                    "thresh": params[5],
                    "recomb_rate": params[4],
                    "mut_rate": mut_rates[m],
                    "spatial_merge_time": params[6],
                    "demo_time": demo_diff,
                    "sim_time": sim_diff,
                    "ts_seed": seed,
                    "mut_seed": mut_seeds[m],
                }, index = [0]
                )
                
                # add demographic summaries to df
                full_df = pd.concat([params_df, sumstats, demo_summary], axis=1)

                # add parameter id
                params_df["param_id"] = f"p{param_id}"
                sumstats["param_id"] = f"p{param_id}"
                demo_summary["param_id"] = f"p{param_id}"

                # add simulation id
                params_df["sim_id"] = f"p{param_id}_{sim_counter}"
                sumstats["sim_id"] = f"p{param_id}_{sim_counter}"
                demo_summary["sim_id"] = f"p{param_id}_{sim_counter}"
    
                param_df_list.append(params_df)
                sumstats_df_list.append(sumstats)
                demo_summary_df_list.append(demo_summary)

                print("dataframes appended")
                sim_counter += 1
                
            
        except Exception as error:
            print(f"Sim with seed {seed} failed. {error}")
    

    # combine parameters into single DataFrame
    param_df_full = pd.concat(param_df_list)
    sumstats_df_full = pd.concat(sumstats_df_list)
    demo_summary_df_full = pd.concat(demo_summary_df_list)

    print("dataframe lists concatenated")

    # write to files
    out_params = f"{out_prefix}_params.csv"
    out_sumstats = f"{out_prefix}_sumstats.csv"
    out_demo = f"{out_prefix}_demo_summary.csv"
    
    if isfile(out_params):
        param_df_full.to_csv(out_params, mode="a", index=False, header=False)
    else:
        param_df_full.to_csv(out_params, index=False)
    
    print("wrote params to file")

    if isfile(out_sumstats):
        sumstats_df_full.to_csv(out_sumstats, mode="a", index=False, header=False)
    else:
        sumstats_df_full.to_csv(out_sumstats, index=False)
    
    print("wrote sumstats to file")

    if isfile(out_demo):
        demo_summary_df_full.to_csv(out_demo, mode="a", index=False, header=False)
    else:
        demo_summary_df_full.to_csv(out_demo, index=False)
    
    print("write demography summaries to file")
    
    print(f"Finished simulation parameter combo id = {param_id}")
    
    # remove unnecessary objects to free up memory and not create any weird situations because tseq and demography are mutable
    del tseq
    del demography
    del demo_summary
    

#########################################
########## Simulations ##################
#########################################

def main():
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

    out_prefix_linear = "sims_2024-02-21_cat_linear_second"
    out_prefix_hinge = "sims_2024-02-18_cat_hinge"

    coords = list(zip(localities.longitude, localities.latitude))

    pop_coords = sampled_cells_to_coords(r2, coords)

    pop_dict_simulation, pop_dict_sim_long, _ = coords_to_sample_dict_empirical(r2, coords)

    max_thresh = max_thresh_from_coords(r2, coords)

    admix_dict = admix_to_deme_dict(admix_rast_path=admix_rast_path, pops_dict=pop_dict_simulation)


    # parameter combinations
    n_p = 10
    
    max_k = np.linspace(100, 5000, n_p,  dtype = int)

    recomb_rate = [1e-9]
    
    mut_rates =  [1e-9]
    
    Nm = np.geomspace(0.1, 5, n_p)

    # spatial merge time
    sp_merge_time_list = [22000]

    thresh = [max_thresh]

    ancestral_size_a1 = np.linspace(5e5, 2e6, n_p, dtype=int)

    ancestral_size_a2 = np.linspace(5e4, 1e6, n_p, dtype=int)

    # since the empirical pi for admix pop 2 (0.110983) is nearly equal to admix pop 3 (0.112128), I'm setting the ancestral size for a3 to be equal to a2

    ancestral_size_list = []

    for i in range(len(ancestral_size_a1)):
        for j in range(len(ancestral_size_a2)):
            if ancestral_size_a1[i] > ancestral_size_a2[j] :
                asize = [ancestral_size_a1[i], ancestral_size_a2[j], ancestral_size_a2[j]]
                ancestral_size_list.append(asize)

    #ancestral_merge_size_list = [[ancestral_size_a1[i]] for i in range(len(ancestral_size_a1))]

    ancestral_merge_time_list = [[1000000]]

    ## missing data percentage to correspond with empirical data
    miss_data = 0.4

    # make list of parameter combinations to iterate over
    lin_sim_params = list(itertools.product(*[
        max_k, 
        Nm, 
        ancestral_size_list, 
       #ancestral_merge_size_list, 
        ancestral_merge_time_list, 
        recomb_rate, 
        [0],
        sp_merge_time_list]))
    # hinge_sim_params = list(itertools.product(*[
    #     max_k, 
    #     Nm, 
    #     ancestral_size_list, 
    #     ancestral_merge_time_list,
    #     recomb_rate, 
    #     thresh,
    #     sp_merge_time_list]))

    # number of simulations per parameter combination
    n_sims_hinge = 5
    # since hinge has an extra parameter that varies (threshold), multiply the number of hinge simulations per parameter combo by the number of threshold parameters explored
    #n_sims_lin = n_thresh * n_sims_hinge
    n_sims_lin = 5

    # iterate over simulations

    # ## I'm separating them into batches because for some reason the Huxley cluster doesn't want to split my jobs across nodes
    mid_ind = len(lin_sim_params) // 2
    lin_sim_params_first = [lin_sim_params[i] for i in range(mid_ind + 1)]
    lin_sim_params_second = [lin_sim_params[i] for i in range(mid_ind, len(lin_sim_params))]

    # continue sims where they left off
    lin_sim_params_second = lin_sim_params_second[4400:]

    # ### linear sims
    num_param_combos = len(lin_sim_params_second)
    sims_per_process = num_param_combos // size
    if num_param_combos % size > 0:
        sims_per_process += 1
    
    first_sim = rank * sims_per_process
    last_sim = first_sim + sims_per_process

    for sim in range(first_sim, last_sim):
        run_sims(params = lin_sim_params_second[sim], 
                 mut_rates = mut_rates,
                 raster = r, 
                 admix_rast = admix_rast, 
                 #merge_time = 22000,
                 pops_dict_sim = pop_dict_simulation, 
                 pops_dict_deme = pop_dict_sim_long, 
                 pops_dict_admix = admix_dict, 
                 demes_coordinates = pop_coords, 
                 n_sims = n_sims_lin, 
                 miss_data = miss_data, 
                 transformation="linear",
                 multi_merge=False,
                 out_prefix = out_prefix_linear)

        
if __name__=="__main__":

    tracemalloc.start()
    start_sims_time = time.time()
    main()
    end_sims_time = time.time()
    total_time_hours = (end_sims_time - start_sims_time) / 3600
    current_mem, max_mem = tracemalloc.get_traced_memory()
    print(f"{current_mem=}, {max_mem=}")
    print(f"Total time = {total_time_hours} hours")
