# Species Distribution Modeling


I modeledthe distributions of the two *Enyalius* study species (*E.
catenatus* and *E. iheringii*) to understand their current distributions
and responses to long-term historical climate change.

I built robust models of their current distributions using CHELSA
bioclims with the objective to project those models to climate time
slices back to 21 kya.

***Note*** The SDMs run using the Maxent binary, which requires Java and
the `rJava` package. Using `rJava` on M-series Apple Macs can be
problematic. I recommend using an older Mac, Windows or Linux machine
for running the SDMs.

``` r
library(sf)
library(spThin)
library(janitor)
library(here)
library(rasterVis)
library(tidyverse)
library(ENMeval)
library(patchwork)
library(terra)
library(usdm)
```

# Current

## Variable selection

I’m choosing variables that encompass the range of the *Enyalius*
species, are relevant to their physiology, and are obtainable at a
reasonable spatial resolution.

Bioclims:  
Using the CHELSA bioclim dataset downloaded from Karger DN, Conrad O,
Böhner J, Kawohl T, Kreft H, Soria-Auza RW, Zimmermann NE, Linder HP,
Kessler M. 2017. Climatologies at high resolution for the earth’s land
surface areas. Scientific Data 4:170122. Downloaded on **2023-11-25**.  
1 km resolution.

For projecting to past climates, I’m using [CHELSA TraCE21k v1.0
Bioclims](https://cp.copernicus.org/preprints/cp-2021-30/cp-2021-30.pdf)
at 100-year intervals back to the LGM. 1 km resolution.

Karger, D.N., Nobis, M.P., Normand, S., Graham, C.H., Zimmermann, N.
(2023) CHELSA-TraCE21k – High resolution (1 km) downscaled transient
temperature and precipitation data since the Last Glacial Maximum.
***Climate of the Past.*** <https://doi.org/10.5194/cp-2021-30>.
Downloaded on **2023-11-25.**

Land Cover:  
Using biome projections from [Costa et
al. 2018](https://onlinelibrary-wiley-com.ezproxy.gc.cuny.edu/doi/full/10.1111/geb.12694)
at 1000-year intervals back to the LGM. 5 km resolution that I’m
downscaling with nearest-neighbor resampling. Downloaded on
**2023-06-08**.

## Variable processing

I’m cropping the variables to a relevant extent before going further
with the individual SDMs. I’m cropping them with a 1 degree buffer
around the entire species group’s range.

First, I’m reading in the localities and plotting them to make sure
there isn’t anything wild. Everything looks reasonable!

``` r
# read in as a data frame first
locs_df <- read_csv(here("analysis", "data", "enyalius_locs.csv")) %>% 
  # the variable names are messy
  janitor::clean_names() %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  # there are some lat/longs with spaces. Need to fix this
  mutate(longitude = str_remove(longitude, "\\s"), 
         latitude = str_remove(latitude, "\\s")) %>% 
  # convert to numeric
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude))

# convert to sf for spatial operations
locs_sf <- st_as_sf(locs_df, 
                    coords = c("longitude", "latitude"),
                    crs = 4326)

# convert to vect for interacting with terra
locs_vec <- vect(locs_sf)

# atlantic forest shapefile for plotting
af <- read_sf(here("analysis", "data", "atlantic_forest", "atlantic_forest.geojson"))

# plot
ggplot() +
  geom_sf(data = st_geometry(af), fill = "gray") +
  geom_sf(data = locs_sf, aes(fill = species), shape = 21) +
  scale_fill_viridis_d() +
  theme_minimal()
```

Creating a convex hull with an 2 degree buffer.

``` r
mcp_all <- st_convex_hull(st_union(locs_sf)) %>%
  st_buffer(dist = units::set_units(2, degree)) %>% 
  terra::vect()

plot(st_geometry(af))
plot(mcp_all, add=TRUE)
plot(locs_vec, add=TRUE)
```

Read in bioclims. Originally I waited to crop, but I folded that into
the initial variable downloads. The rasters going into
`cropped_predictors` are the same as I’m reading in, just with different
names. They’re also being written as a single file, rather than multiple
files.

``` r
bioclims <- terra::rast(list.files(here("analysis", "data", "current_climate_chelsa"), full.names = TRUE))

# get the only the bioclim labels
names(bioclims) <- str_extract(names(bioclims), "bio[\\d+]*")

plot(bioclims$bio1)
plot(st_geometry(af), add=TRUE)
plot(mcp_all, add=TRUE)
plot(locs_vec, add=TRUE)
```

Write the cropped climate layers to file for use in the SDMs.

``` r
terra::writeRaster(bioclims, here("analysis", "output", "cropped_predictors", "bioclims.tif"))
```

## General SDM steps

For each SDM, I’m going to perform the following steps:

- spatially thin localities
  - 20 km buffer to reduce spatial autocorrelation  
- crop the new predictors according to a minimum convex hull of the
  localities for each species with a 0.5 degree buffer
  - I tried point-buffers, but they resulted in worse models  
- extract background environmental data (10,000 points) for predictor
  correlation and modeling  
- correlate predictors  
- remove predictors w/ r \> 0.75, prioritizing ecologically relevant
  variables  
- Use Maxent for modeling
  - LQH feature classes to keep models relatively simple  
  - regularization multipliers from 0.5 to 5.0 in 0.5 increments to test
    a wide range of regularization.
  - leave-one-out cross validation for model selection due to a low
    number of localities  
  - select models first by AICc (model fit), followed by ommission error
    rate (prediction)

In addition, I’m reading/writing data for each species in isolation, so
the code for one species isn’t dependent on that for another, or for
what I did during variable processing.

### E. catenatus

Note: putative “catenatus 2” species from Mariana’s phylogenetic work
have been removed.

Reading in data for plotting.

``` r
# atlantic forest shapefile
af <- read_sf(here("analysis", "data", "atlantic_forest", "atlantic_forest.geojson"))
```

#### Spatial thin

Read and filter localities.

``` r
locs_cat <- read_csv(here("analysis", "data", "enyalius_locs.csv")) %>% 
  # the variable names are messy
  janitor::clean_names() %>%
  filter(species == "catenatus",
    !is.na(latitude), !is.na(longitude),
    # remove duplicates
    !duplicated(latitude)) %>%
  # there are some lat/longs with spaces. Need to fix this
  mutate(longitude = str_remove(longitude, "\\s"), 
         latitude = str_remove(latitude, "\\s")) %>% 
  # convert to numeric
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>% 
  # one locality falls into the ocean when using the aggregated data. Have to scooch it over by 0.05 degrees to overlap with the correct cell
  mutate(longitude = if_else(longitude == -39.06, longitude - 0.05, longitude)) %>% 
  # convert to sf for spatial operations
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326,
           remove = FALSE)

# plot localities
plot(st_geometry(af))
plot(st_geometry(locs_cat), add = TRUE)
```

Spatial thin. I’m using a 20 km buffer to reduce the impact of spatial
autocorrelation.

``` r
set.seed(394833)

#run spthin algorithm. This returns 100 possible combinations of removed localities
output <-
  spThin::thin(
    locs_cat,
    'latitude',
    'longitude',
    'species',
    thin.par = 20,
    reps = 100,
    locs.thinned.list.return = TRUE,
    write.files = FALSE,
    verbose = FALSE
  )

# I want to maximize the # of localities returned  
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))

# if there are multiple iterations with the max # of localities, pick one
maxThin <-
  output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]

# subset locs to match only thinned locs
locs_cat <- locs_cat[as.numeric(rownames(maxThin)), ]

# get the unused locs as a testing data set
# this probably isn't useful since they will overlap heavily with the training data, but seems like another piece of info to look at
test_cat <- locs_cat[-as.numeric(rownames(maxThin)), ]

plot(st_geometry(af))
plot(st_geometry(locs_cat), add=TRUE)
```

Write thinned localities to file

``` r
# Write to file
st_write(locs_cat, here("analysis", "output", "thinned_localities", "catenatus_thinned.gpkg"))
```

#### Crop environment

First, I need to read in the environmental data.

``` r
bioclims <- terra::rast(here("analysis", "output", "cropped_predictors", "bioclims.tif"))
```

Cropping and combining variables for analysis.

``` r
# for projection
mcp_cat <- st_convex_hull(st_union(locs_cat)) %>%
  st_buffer(dist = units::set_units(0.5, degree)) %>% 
  terra::vect()

bioclims_cat <- bioclims %>%
  terra::crop(mcp_cat) %>% 
  terra::mask(mcp_cat)

predictors_cat <- c(bioclims_cat)

plot(bioclims_cat[[1]])
plot(locs_cat, add=TRUE)
```

#### Predictor correlation

Sample 10000 background points. Only returned 2355 points.

``` r
set.seed(1988888)

# for variable correlations
bg_envt_cat <- terra::spatSample(predictors_cat, 10000, 
                               warn=TRUE, 
                               na.rm = TRUE, 
                               as.df = TRUE,
                               xy = TRUE)

# for use in ENMevaluate
bg_coords_cat <- bg_envt_cat[,c("x", "y")]
```

Next, I’ll extract the values for the background points and perform
variance inflation factor stepwise selection with a VIF threshold of 10.

``` r
# extract values
bg_corr_cat <- bg_envt_cat %>% select(-x, -y)

usdm::vifstep(bg_corr_cat, th=10)
```

The final variable list: BIO3, BIO4, BIO8, BIO13, BIO15, BIO18

``` r
# label: pred-cat

predictors_cat <- predictors_cat[[c("bio3", "bio4", "bio8", "bio13", "bio15", "bio18")]]
```

#### Maxent model

I’m using a jackknifing model evaluation approach since I only have 26
observations and a previous attempt for spatial CV led to wonky
evaluation metrics. Additionally, spatial CV can lead to worse
predictions than non-spatial CV: [Wadoux et
al. 2021](https://www-sciencedirect-com.ezproxy.gc.cuny.edu/science/article/pii/S0304380021002489).

``` r
set.seed(7990777)
coords_cat <- st_coordinates(locs_cat)
colnames(coords_cat) <- c("x", "y")

folds_cat <- ENMeval::get.jackknife(occ = coords_cat, 
                            bg = bg_coords_cat)
```

Run the model. Predictions are clamped to prevent extrapolation.

``` r
set.seed(34622)

# the vector of regularization multipliers to test
rms <- seq(0.5, 5, 0.5)

# convert the terra object to a raster stack for use in EMNeval
predictors_cat_stack <- raster::stack(predictors_cat)

# iterate model building over all chosen parameter settings
sdm_cat <-
  ENMeval::ENMevaluate(
    occs = coords_cat,
    envs = predictors_cat_stack,
    bg.coords = bg_coords_cat,
    RMvalues = rms,
    fc = c('L', 'LQ'),
    # clamping to prevent model extrapolation
    doClamp = TRUE,
    taxon.name = "catenatus",
    partitions = "user",
    user.grp = folds_cat,
    bg.grp,
    clamp = TRUE,
    algorithm = "maxent.jar",
    parallel = TRUE,
    numCores = 6
  )
```

``` r
# write the model to file
write_rds(sdm_cat, here("analysis", "output", "sdm_models", "sdm_catenatus.rds"))
```

##### Model evaluation

Let’s take a look at the model results.

``` r
eval_table_cat <- sdm_cat@results
eval_mods_cat <- sdm_cat@models

names(eval_mods_cat) <-
  str_replace_all(names(eval_mods_cat), "\\.", "\\_")
```

Select the final model. First I’m looking at plots of model evaluation
stats to get an idea of the spread of stats.

``` r
daic_cat <- ggplot(data = eval_table_cat, aes(x = rm, y = delta.AICc, color = fc)) +
  geom_point() +
  scale_color_viridis_d() +
  theme_bw()

or_cat <- ggplot(data = eval_table_cat, aes(x = rm, y = or.10p.avg, color = fc)) +
  geom_point() +
  scale_color_viridis_d() +
  theme_bw()

dauc_cat <- ggplot(data = eval_table_cat, aes(x = rm, y = auc.diff.avg, color = fc)) +
  geom_point() +
  scale_color_viridis_d() +
  theme_bw()

(daic_cat + or_cat) / (dauc_cat)
```

Now I’m going to take a look at tables of delta AICc, omission rate, and
AUC to see how close the models are. The model with the lowest AICc has
a regularization multiplier of 1 and L feature class.

``` r
eval_table_cat %>% 
  select(delta.AICc, AICc, or.10p.avg, or.mtp.avg, auc.diff.avg, auc.val.avg, rm, fc) %>%
  arrange(delta.AICc) %>% 
  head(10) %>% 
  knitr::kable()
```

Select model.

``` r
mod_cat <- eval_mods_cat$rm_0_5_fc_L
opt_seq_tune_cat <- "rm.0.5_fc.L"
```

Plot the variable contributions and response curves. The most important
variables are bio15, bio18, and bio13. All precipitation!

``` r
plot(mod_cat)
```

``` r
png("variable_contribution_catenatus.png")
plot(mod_cat)
dev.off()

file.copy(here("variable_contribution_catenatus.png"), here("analysis", "output", "sdm_response_curves_variable_importance", "variable_contribution_catenatus.png"))
file.remove(here("variable_contribution_catenatus.png"))
```

Suitability increases with precipitation amount of the wettest month,
decreases with precipitation seasonality, and increases with mean
monthly precipitation of the warmest quarter.

``` r
dismo::response(mod_cat)
```

``` r
png("response_curves_catenatus.png")
dismo::response(mod_cat)
dev.off()

file.copy(here("response_curves_catenatus.png"), here("analysis", "output", "sdm_response_curves_variable_importance", "response_curves_catenatus.png"))
file.remove(here("response_curves_catenatus.png"))
```

##### Project

I’m projecting the model to the study area extent.

``` r
pred_cat <- ENMeval::eval.predictions(sdm_cat)[[opt_seq_tune_cat]]
plot(pred_cat)
plot(st_geometry(af), add = TRUE)
plot(st_geometry(locs_cat), pch = 21, bg = alpha("lightgray", 0.5), add = TRUE)
```

``` r
pred_cat_thresh <- pred_cat$rm.0.5_fc.L
# observed suitabilities
obs_suit_cat <- terra::extract(pred_cat$rm.0.5_fc.L, locs_cat)
# minimum training presence
min_suit_cat <- min(obs_suit_cat)

pred_cat_thresh[pred_cat_thresh >= min_suit_cat] <- 1
pred_cat_thresh[pred_cat_thresh < min_suit_cat] <- 0

plot(pred_cat_thresh)
plot(st_geometry(af), add = TRUE)
plot(st_geometry(locs_cat), pch = 21, bg = alpha("lightgray", 0.5), add = TRUE)
```

### E. iheringii

Reading in data for plotting.

``` r
# atlantic forest shapefile
af <- read_sf(here("analysis", "data", "atlantic_forest", "atlantic_forest.geojson"))
```

#### Spatial thin

Read and filter localities.

``` r
locs_ihe <- read_csv(here("analysis", "data", "enyalius_locs.csv")) %>% 
  # the variable names are messy
  janitor::clean_names() %>%
  filter(species == "iheringii",
    !is.na(latitude), !is.na(longitude),
    # remove duplicates
    !duplicated(latitude)
    ) %>%
  # there are some lat/longs with spaces. Need to fix this
  mutate(longitude = str_remove(longitude, "\\s"), 
         latitude = str_remove(latitude, "\\s")) %>% 
  # convert to numeric
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>% 
  # some localities fall into the ocean when the rasters are aggregated. They need to be bumped by 0.05 degrees to overlap their respective environmental layer. Given the coarseness of the analysis, this does not bias results. 
  mutate(
    longitude = case_when(
      longitude == -47.01956 ~ longitude - 0.05,
      longitude == -47.93278 ~ longitude - 0.05,
      longitude == -47.93139 ~ longitude - 0.05,
      longitude == -48.03889 ~ longitude - 0.05,
      .default = longitude
    ),
    latitude = case_when(latitude == -23.595 ~ latitude + 0.05,
                         .default = latitude)
  ) %>% 
  # convert to sf for spatial operations
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326,
           remove = FALSE)

# plot localities
plot(st_geometry(af))
plot(st_geometry(locs_ihe), add = TRUE)
```

Spatial thin. I’m using a 20 km buffer

``` r
set.seed(2333)

#run spthin algorithm. This returns 100 possible combinations of removed localities
output <-
  spThin::thin(
    locs_ihe,
    'latitude',
    'longitude',
    'species',
    thin.par = 20,
    reps = 100,
    locs.thinned.list.return = TRUE,
    write.files = FALSE,
    verbose = FALSE
  )

# I want to maximize the # of localities returned  
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))

# if there are multiple iterations with the max # of localities, pick one
maxThin <-
  output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]

# subset locs to match only thinned locs
locs_ihe <- locs_ihe[as.numeric(rownames(maxThin)), ]

# get the unused locs as a testing data set
# this probably isn't useful since they will overlap heavily with the training data, but seems like another piece of info to look at
test_ihe <- locs_ihe[-as.numeric(rownames(maxThin)), ]

plot(st_geometry(af))
plot(st_geometry(locs_ihe), add=TRUE)
```

Write thinned localities to file

``` r
# Write to file
st_write(locs_ihe, here("analysis", "output", "thinned_localities", "iheringii_thinned.gpkg"),
         delete_dsn = TRUE)
```

#### Crop environment

First, I need to read in the environmental data.

``` r
bioclims <- terra::rast(here("analysis", "output", "cropped_predictors", "bioclims.tif"))
```

Cropping and combining variables for analysis. Cropping to 0.5 degree
buffer to accomodate reasonable dispersal

``` r
mcp_ihe <- st_convex_hull(st_union(locs_ihe)) %>%
  st_buffer(dist = units::set_units(0.5, degree)) %>% 
  terra::vect()

bioclims_ihe <- bioclims %>%
  terra::crop(mcp_ihe) %>% 
  terra::mask(mcp_ihe)

plot(bioclims_ihe[[1]])
plot(st_geometry(locs_ihe), add = TRUE)
plot(st_geometry(locs_ihe), add = TRUE)

predictors_ihe <- c(bioclims_ihe)
```

#### Predictor correlation

Sample 10000 background points. Only 6515 were able to be sampled.

``` r
set.seed(7488)

# for variable correlations
bg_envt_ihe <- terra::spatSample(predictors_ihe, 10000, 
                               warn=TRUE, 
                               na.rm = TRUE, 
                               as.df = TRUE,
                               xy = TRUE)

# for use in ENMevaluate
bg_coords_ihe <- bg_envt_ihe[,c("x", "y")]
```

Next, I’ll extract the values for the background points and perform
variance inflation factor stepwise selection with a VIF threshold of 10.

``` r
# extract values
bg_corr_ihe <- bg_envt_ihe %>% select(-x, -y)

usdm::vifstep(bg_corr_ihe, th=10)
```

The final variable list: BIO3, BIO7, BIO8, BIO9, BIO14, BIO18

``` r
predictors_ihe <- predictors_ihe[[c("bio3", "bio7", "bio8", "bio9", "bio14", "bio18")]]
```

#### Maxent model

I’m using a jackknife model evaluation approach since I only have 26
observations.

``` r
set.seed(3398737)
coords_ihe <- st_coordinates(locs_ihe)
colnames(coords_ihe) <- c("x", "y")

folds_ihe <- ENMeval::get.jackknife(occ = coords_ihe, 
                            bg = bg_coords_ihe)
```

Run the model. Predictions are clamped to prevent extrapolation.

``` r
set.seed(19999923)

# the vector of regularization multipliers to test
rms <- seq(0.5, 5, 0.5)

# convert the terra object to a raster stack for use in EMNeval
predictors_ihe_stack <- raster::stack(predictors_ihe)

# iterate model building over all chosen parameter settings
sdm_ihe <-
  ENMeval::ENMevaluate(
    occs = coords_ihe,
    envs = predictors_ihe_stack,
    bg.coords = bg_coords_ihe,
    RMvalues = rms,
    fc = c('L', 'LQ'),
    # clamping to prevent model extrapolation
    doClamp = TRUE,
    taxon.name = "iheringii",
    partitions = "user",
    # going to do a separate run with the final model on testing data
    #occs.testing = test_ihe,
    user.grp = folds_ihe,
    clamp = TRUE,
    algorithm = "maxent.jar",
    parallel = TRUE,
    numCores = 6
  )
```

``` r
# write the model to file
write_rds(sdm_ihe, here("analysis", "output", "sdm_models", "sdm_iheringii.rds"))
```

##### Model evaluation

Let’s take a look at the model results.

``` r
eval_table_ihe <- sdm_ihe@results
eval_mods_ihe <- sdm_ihe@models

names(eval_mods_ihe) <-
  str_replace_all(names(eval_mods_ihe), "\\.", "\\_")
```

Select the final model. First I’m looking at plots of model evaluation
stats to get an idea of the spread of stats.

``` r
daic_ihe <- ggplot(data = eval_table_ihe, aes(x = rm, y = delta.AICc, color = fc)) +
  geom_point() +
  scale_color_viridis_d() +
  theme_bw()

or_ihe <- ggplot(data = eval_table_ihe, aes(x = rm, y = or.10p.avg, color = fc)) +
  geom_point() +
  scale_color_viridis_d() +
  theme_bw()

dauc_ihe <- ggplot(data = eval_table_ihe, aes(x = rm, y = auc.diff.avg, color = fc)) +
  geom_point() +
  scale_color_viridis_d() +
  theme_bw()

(daic_ihe + or_ihe) / (dauc_ihe)
```

Now I’m going to take a look at tables of delta AICc, omission rate, and
AUC to see how close the models are. All of the top models have the same
omission rate, similar AUC, and are either LQ or L. I chose the simplest
model within 2 AICc of the top model, which was the rm 2.5 and linear
feature class model. The top was a rm 2.5 linear quadratic, but even
marginally simpler is better within reasonable limits. After skimming
through them, they are all very similar models.

``` r
eval_table_ihe %>% 
  select(delta.AICc, AICc, or.10p.avg, or.mtp.avg, auc.diff.avg, auc.val.avg, rm, fc) %>%
  arrange(delta.AICc) %>% 
  head(10) %>% 
  knitr::kable()
```

Select model.

``` r
mod_ihe <- eval_mods_ihe$rm_2_fc_LQ
opt_seq_tune_ihe <- eval_table_ihe$tune.args[eval_table_ihe$tune.args == "rm.2_fc.LQ"]
```

Variable importance.

``` r
plot(mod_ihe)
```

``` r
png("variable_contribution_iheringii.png")
plot(mod_ihe)
dev.off()

file.copy(here("variable_contribution_iheringii.png"), here("analysis", "output", "sdm_response_curves_variable_importance", "variable_contribution_iheringii.png"))
file.remove(here("variable_contribution_iheringii.png"))
```

Plot the response curves. In order: bio3, bio7^2, bio8, bio9, bio14^2,
bio18.

``` r
dismo::response(mod_ihe)
```

``` r
png("response_curves_iheringii.png")
dismo::response(mod_ihe)
dev.off()

file.copy(here("response_curves_iheringii.png"), here("analysis", "output", "sdm_response_curves_variable_importance", "response_curves_iheringii.png"))
file.remove(here("response_curves_iheringii.png"))
```

##### Project

I’m projecting the model to the study area extent. It’s not predicting
strong suitability for the southern localities. I’ll have to look at the
genetic structure and locality info more closely to see what’s up with
them.

``` r
pred_ihe <- ENMeval::eval.predictions(sdm_ihe)[[opt_seq_tune_ihe]]
plot(pred_ihe)
plot(st_geometry(af), add = TRUE)
plot(st_geometry(locs_ihe), pch = 21, bg = alpha("lightgray", 0.5), add = TRUE)
```

``` r
pred_ihe_thresh <- pred_ihe$rm.2_fc.LQ
# observed suitabilities
obs_suit_ihe <- terra::extract(pred_ihe$rm.2_fc.LQ, locs_ihe)
# minimum training presence
min_suit_ihe <- min(obs_suit_ihe)

pred_ihe_thresh[pred_ihe_thresh >= min_suit_ihe] <- 1
pred_ihe_thresh[pred_ihe_thresh < min_suit_ihe] <- 0

plot(pred_ihe_thresh)
plot(st_geometry(af), add = TRUE)
plot(st_geometry(locs_ihe), pch = 21, bg = alpha("lightgray", 0.5), add = TRUE)
```
