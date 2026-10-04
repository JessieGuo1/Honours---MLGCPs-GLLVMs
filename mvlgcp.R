# Get the job array
load("job_array.RDATA")

# determine job number from pbs script
job = 1

## MODEL SETTINGS ##############################################################

# set the number of latent factor we want to use
d = 1
# set approximate number of basis function nodes to use
basis.dim = 43
# get the predictor set to be modelled
# pred_set <- as.character(tab$p[tab$job == job])
# get the basis function type to be used
# bf.type <- as.character(tab$bf.type[tab$job == job])
################################################################################

#DIFFERENT DATASET
library(disdat)
library(mgcv)
library(sf)
library(spatstat)
library(fields)
library(tmap)
library(grid)
library(units)
library(stringr)
library(dplyr)
library(glmnet)
library(scampr)
library(glmmTMB)
set.seed(0)

#load data
po = disPo("NSW")
bg = disBg("NSW")
border_lonlat = disBorder("NSW")

#convert presence/background data frames to sf objects using their lon/lat coords (EPSG:4326)
po_lonlat = st_as_sf(po, coords = c("x", "y"), crs = 4326, remove = FALSE) #remove = FALSE keeps x, y
bg_lonlat = st_as_sf(bg, coords = c("x", "y"), crs = 4326, remove = FALSE)

target_crs = 3577 # Australian Albers, coordinates in metres

#project border, presence, and background points from lon/lat into Albers metres
border_proj = st_transform(border_lonlat, crs = target_crs)
po_proj = st_transform(po_lonlat, crs = target_crs) 
bg_proj = st_transform(bg_lonlat, crs = target_crs)

#extract projected coordinates as a plain matrix
coords_po = st_coordinates(po_proj)
coords_bg = st_coordinates(bg_proj)

#drop geometry column --> plain data frame
po_model = st_drop_geometry(po_proj)
bg_model = st_drop_geometry(bg_proj)

#rescale projected coordinates from metres to kilometres, for numerical stability
po_model$x = coords_po[, "X"]/1000
po_model$y = coords_po[, "Y"]/1000

bg_model$x <- coords_bg[, "X"]/1000
bg_model$y <- coords_bg[, "Y"]/1000

#covariates of data
covar_names = c("cti", "disturb", "mi", "rainann", "raindq", "rugged",
                "soildepth", "soilfert", "solrad", "tempann", "tempmin", "topo", "vegsys")

#keep rows with complete covariate data
po_keep = complete.cases(po_model[, covar_names])
bg_keep = complete.cases(bg_model[, covar_names])

po_lonlat = po_lonlat[po_keep, ]
bg_lonlat = bg_lonlat[bg_keep, ]

#total study region area (km^2), used to compute quadrature weights
domain_area_km2 <- as.numeric(set_units(sum(st_area(border_proj)), km^2))

load("bg_model.RDATA")
quad_wt <- domain_area_km2 / nrow(bg_model) #each background/quadrature point represents domain_area_km2

#processed dataframes
presence_df = data.frame(
  spid = po_model$spid, 
  x = po_model$x, y = po_model$y,
  po_model[, covar_names],
  presence = 1, 
  quad.size = 0) #quadrature weights not used at the point events

quad = data.frame(
  x = bg_model$x, y = bg_model$y,
  bg_model[, covar_names],
  presence = 0,
  quad.size = quad_wt) #area represented per quad point

species = unique(po_model$spid)
m = length(species) #no. of species
q <- nrow(quad) #no. of quadrature points

#replicate background points once per species for multi-species scheme)
quad <- quad[rep(seq_len(q), m), ]
quad$spid <- rep(species, each = q) #label each replicate block with its species id

#combine presence, quadrature rows into long-format dataframe
dat <- rbind(presence_df, quad)
n <- nrow(dat)

#pt/wt - standard DWPR (down-weighted Poisson regression) construction
#pt = 1 - presence rows, 0 - quadrature rows
#wt = tiny fixed weight for presence rows, area-based weight for quadrature rows
dat$pt = dat$presence
dat$pt = dat$presence
dat$wt = ifelse(dat$presence == 1, 1e-6, dat$quad.size)

cont_vars <- setdiff(covar_names, "vegsys")
dat_scaled = scale(dat[, cont_vars])
dat[, cont_vars] = as.data.frame(dat_scaled) 

#evaluate basis functions at every data point
bfs_elliot <- make_basis(k = basis.dim, quad, from.package = "scampr")

# need to prune in the case of scampr bfs
if (is(bfs_elliot, "scampr.bf")) {
  tmp <- bf_matrix(bfs_elliot, rbind(presence_df[,c("x", "y")], quad[,c("x", "y")])) # via the data?
  # tmp <- bf_matrix(bfs, pres[,c("xm", "ym")]) # via the data?
  bf.check <- apply(tmp, 2, function(x){sum(x > 0)})
  bfs_elliot <- bfs_elliot[bf.check != 0, ]
  rm(tmp, bf.check)
}

dat_rename = dat %>% 
  mutate(spid = as.numeric(str_extract(spid, "\\d+")))
#  mutate(Zu = sweep(Z_star, 2, u, "*"))#
#  

attach(dat_rename)

chosen_spid = unique(spid)
n_species = length(chosen_spid)

#rm(pres, quad)
#gc()

time_i = system.time({try(assign("mB1", mvlgcp(pt ~ 1 + tempann + rainann + rugged,
                         data = dat_rename, weights = wt, basis.functions = bfs_elliot, response.id = spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))})

# save if the model fitted
if (exists("mB1")) {
  print("mB1 converged without start parameters (attempt 1)")
#  mB1 <- rm.env(mB1)
  mB1$start.pars <- "cold"
  save(mB1, file = paste0(getwd(), "/models/mB1_", job, ".RDATA"))
}

pred_mlgcp = predict(mB1, type = "response")

loglik_new_mlgcp = sum(dat_rename$pt * log(pred_mlgcp) - dat_rename$wt * pred_mlgcp)
save(loglik_new_mlgcp, file = "loglik_new_mlgcp.RDATA")

time_mlgcp = as.numeric(time_i[3])
time_mlgcp
save(time_mlgcp, file = "time_mlgcp.RDATA")


#prediction grid coordinates + covariates
load("grid_xy.RDATA")

#fixed-effect design matrix at grid locations
X_grid = model.matrix( ~tempann + rainann + rugged, data = grid_xy) #18277 x 4

#evaluate basis functions at grid locations, using the same bfs as model fitting
if (is(bfs_elliot, "scampr.bf")) {
  Z_grid_star <- bf_matrix(bfs_elliot, grid_xy[, c("x", "y")])
}

#spatial term at grid locations, scaled by the fitted shared u
#Zu_grid = sweep(Z_grid_star, 2, u, "*")
pars = mB1$obj$env$parList()

B = matrix(pars$b, nrow = ncol(Z_grid_star), ncol = length(chosen_spid))
ZB_grid = Z_grid_star%*%B

target_species_id = 9
species_idx = which(chosen_spid == target_species_id)

beta = pars$beta
#predicted intensity for every species at every grid cell, then select just one species' column
z_matrix = exp(X_grid%*%beta + ZB_grid)

grid_xy$z = z_matrix[, species_idx]

#builds intensity map - study region border + fitted intensity surface + this species' presence points
load("m_intensity_map.RDATA")

#spatial polygons matching grid_xy locations
load("grid_cells_proj.RDATA") 
#attach predicted intensity to grid cell geometries, transform to lon/lat for plotting
grid_pred_proj = st_sf(intensity = grid_xy$z,
                       log_intensity = log(grid_xy$z),
                       geometry = grid_cells_proj)
grid_pred_lonlat = st_transform(grid_pred_proj, target = 4326)

#pull this species' own presence points, convert km back to metres, then to lon/lat for plotting
species_points = dat_rename %>% filter(spid == target_species_id, pt == 1)
species_points_metres = species_points
species_points_metres$x = species_points$x * 1000
species_points_metres$y = species_points$y * 1000

species_points_sf = st_as_sf(species_points_metres, coords = c("x", "y"), crs = target_crs)
species_points_lonlat = st_transform(species_points_sf, crs = 4326)

glmnet_mlgcp_map = m_intensity_map(grid_pred_lonlat, species_points_lonlat, species_label = 9)
print(glmnet_mlgcp_map)
save(glmnet_mlgcp_map, file = "glmnet_mlgcp_map.RDATA")










# with dual fields #############################################################

# first try to fit with cold starts for parameters
try(assign("mB2", mvlgcp(occ ~ (1 | spid) + mi + tempann + (0 + mi | spid) + (0 + tempann | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both")
))
if (exists("mB2")) {
  print(paste0("mB2 converged without start parameters (attempt 1)", if(is.na(logLik(mB2))) {" however, the likelihood was NA"} ))
  mB2$start.pars <- "cold"
  if (is.na(logLik(mB2))) {
    rm(mB2)
  }
}

# next try with mB1 as starts for parameters
if (exists("mB1") & !exists("mB2")) {
  start.pars <- lapply(split(mB1$fit$parfull, names(mB1$fit$parfull)), unname)
  start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * m))
  start.pars$theta <- c(start.pars$theta, rep(-3, m))
  try(assign("mB2", mvlgcp(occ ~ (1 | spid) + mi + tempann + (0 + mi | spid) + (0 + tempann | spid),
                           data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                           control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both",
                           start = start.pars)
  ))
  rm(start.pars)
  # check that the model has a log-likelihood
  if (exists("mB2")) {
    print(paste0("mB2 fitted with starting parameters from mB1 (attempt 2)", if(is.na(logLik(mB2))) {" however, the likelihood was NA"} ))
    mB2$start.pars <- "warm"
    # if (is.na(logLik(mB2))) {
    #   rm(mB2)
    # }
  }
}

# save if the model fitted
if (exists("mB2")) {
  mB2 <- rm.env(mB2)
  save(mB2, file = paste0(getwd(), "/models/mB2_", job, ".RDATA"))
}

################################################################################
# Fit the model including quadratic env. terms #################################
################################################################################

# with correlated fields only ##################################################

try(assign("mC1", mvlgcp(occ ~ (1 | spid) + mi + tempann + mi2 + tempann2 + (0 + mi |  spid) +
                           (0 + tempann | spid) + (0 + mi2 | spid) + (0 + tempann2 | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))

# save if the model fitted
if (exists("mC1")) {
  mC1 <- rm.env(mC1)
  mC1$start.pars <- "cold"
  print("mC1 converged without start parameters (attempt 1)")
  save(mC1, file = paste0(getwd(), "/models/mC1_", job, ".RDATA"))
}

# with dual fields #############################################################

# first try to fit with cold starts for parameters
try(assign("mC2", mvlgcp(occ ~ (1 | spid) + mi + tempann + mi2 + tempann2 + (0 + mi |  spid) +
                           (0 + tempann | spid) + (0 + mi2 | spid) + (0 + tempann2 | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both")
))
if (exists("mC2")) {
  print(paste0("mC2 converged without start parameters (attempt 1)", if(is.na(logLik(mC2))) {" however, the likelihood was NA"} ))
  mC2$start.pars <- "cold"
  if (is.na(logLik(mC2))) {
    rm(mC2)
  }
}

# next try with mC1 as starts for parameters
if (exists("mC1") & !exists("mC1")) {
  start.pars <- lapply(split(mC1$fit$parfull, names(mC1$fit$parfull)), unname)
  start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * m))
  start.pars$theta <- c(start.pars$theta, rep(-3, m))
  try(assign("mC2", mvlgcp(occ ~ (1 | spid) + mi + tempann + mi2 + tempann2 + (0 + mi |  spid) +
                             (0 + tempann | spid) + (0 + mi2 | spid) + (0 + tempann2 | spid),
                           data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                           control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both",
                           start = start.pars)
  ))
  rm(start.pars)
  # check that the model has a log-likelihood
  if (exists("mC2")) {
    print(paste0("mC2 fitted with starting parameters from mC1 (attempt 2)", if(is.na(logLik(mC2))) {" however, the likelihood was NA"} ))
    mC2$start.pars <- "warm"
    if (is.na(logLik(mC2))) {
      rm(mC2)
    }
  }
}

# save if the model fitted
if (exists("mC2")) {
  mC2 <- rm.env(mC2)
  save(mC2, file = paste0(getwd(), "/models/mC2_", job, ".RDATA"))
}

