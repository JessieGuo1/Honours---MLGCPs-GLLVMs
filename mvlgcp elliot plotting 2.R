# set approximate number of basis function nodes to use
basis.dim = 50

library(glmmTMB)

# get the data
source("get_and_format_flora_data.R")
rm(dat_pa)

# get a vector of the factor terms
facts <- c("disturb", "soilfert", "vegsis")

# load the species information
load("flora_species.RDATA")

# get the groups to which the flora species belongs to subset the data
flora_groups <- unique(sp.info$group)

# subset the presence-only data (in the case we are using restricted species)
pres <- pres[pres$group %in% flora_groups, ]

# want to use species with at least 10 presences
pres <- pres[pres$spid %in% names(table(pres$spid))[table(pres$spid) >= 10], ]

# get the species
species <- unique(pres$spid)

# assign some numbers
m <- length(species) # responses
q <- nrow(quad) # quadrature points

# add in the quadrature for each species (replicates of the supplied background points)
tmp <- quad[rep(seq_len(q), m), ]
tmp$spid <- rep(species, each = q) # assign the species id's
dat <- rbind(pres, tmp)
n <- nrow(dat)
rm(tmp)

# set the coordinates to x and y so we can just use default labels
dat$lon <- dat$x; dat$lat <- dat$y; dat$x <- dat$xm; dat$y <- dat$ym; dat$xm <- NULL; dat$ym <- NULL
domain$lon <- domain$x; domain$lat <- domain$y; domain$x <- domain$xm; domain$y <- domain$ym; domain$xm <- NULL; domain$ym <- NULL

# set the basis functions

bfs <- make_basis(k = basis.dim, domain, from.package = "scampr")

# need to prune in the case of scampr bfs
if (is(bfs, "scampr.bf")) {
  tmp <- bf_matrix(bfs, rbind(pres[,c("xm", "ym")], quad[,c("xm", "ym")])) # via the data?
  # tmp <- bf_matrix(bfs, pres[,c("xm", "ym")]) # via the data?
  bf.check <- apply(tmp, 2, function(x){sum(x > 0)})
  bfs <- bfs[bf.check != 0, ]
  rm(tmp, bf.check)
}

gc()


################################################################################
# Fit the model including linear env terms #####################################
################################################################################
set.seed(0)
# with correlated fields only ##################################################
time_i = system.time({try(assign("mB2", mvlgcp(occ ~ 0 + spid + spid:mi + spid:tempann,
                                               data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = 2,
                                               control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))})

# save if the model fitted
if (exists("mB2")) {
  print("mB2 converged without start parameters (attempt 1)")
  #  mB2 <- rm.env(mB2)
  mB2$start.pars <- "cold"
  save(mB2, file = paste0(getwd(), "/models/mB2.RDATA"))
}


pred_mlgcp = predict(mB2, type = "response", re.form = NULL)

loglik_new_mlgcp_multi_u = sum(dat$occ * log(pred_mlgcp) - dat$wt * pred_mlgcp)
save(loglik_new_mlgcp_multi_u, file = "results/loglik_new_mlgcp_multi_u.RDATA")

time_mlgcp_multi_u = as.numeric(time_i[3])
time_mlgcp_multi_u
save(time_mlgcp_multi_u, file = "results/time_mlgcp_multi_u.RDATA")

#Using tmap
target_species_id <- "nsw19"
Lambda = get_loadings(mB2)
species_idx = which(species == target_species_id)

get_field_fixed = function(object, newdata, which.response = NULL, combine.response.fields = TRUE) {
  if (is.null(object$basis.functions)) stop("Supplied fitted model does not include basis functions")
  
  if (missing(newdata)) {
    warning("'newdata' is not provided. Fitted values will be used...")
    newdata = object$frame$coords
  }
  
  if (ncol(newdata) > 2) {
    if (!all(attr(object$basis.functions, "coord.names") %in% colnames(newdata))) {
      stop("'coord.names' associated with the fitted model's basis functions are not found within the newdata supplied")
    } else {
      newdata = newdata[, attr(object$basis.functions, "coord.names")]
    }
  }
  
  X = bf_matrix(object$basis.functions, newdata)
  b = get_bfs_coeff(object)
  fld = X %*% b
  
  if (object$n_responses == 1) {
    Lambda = NULL
  } else {
    Lambda = get_loadings(object)
  }
  
  if (!is.null(which.response)) {
    if (object$n_responses == 1) {
      warning("'which.response' is not applicable to a univariate model... returning the single latent field")
      Lambda = NULL
    } else if (object$n_responses < length(which.response)) {
      stop("more fields requested through 'which.response' than responses in the fitted model. Please check and retry")
    } else {
      
      if (object$mv.flds %in% c("correlated", "both")) {
        fld_cor = fld[, grepl("d", colnames(fld), fixed = TRUE), drop = FALSE] %*% t(Lambda[which.response, , drop = FALSE])
      } else {
        fld_cor = NULL
      }
      
      if (!is.null(fld_cor)) colnames(fld_cor) = paste0("m", which.response)
      
      if (object$mv.flds %in% c("independent", "both")) {
        if (length(which.response) > 1) {
          fld_ind = fld[, paste0("m", which.response), drop = FALSE]
        } else {
          fld_ind = as.matrix(fld[, paste0("m", which.response)])
        }
      } else {
        fld_ind = NULL
      }
      
      if (!is.null(fld_ind)) colnames(fld_ind) = paste0("m", which.response, "_ind")
      
      if (combine.response.fields && !is.null(fld_cor) && !is.null(fld_ind)) {
        fld = fld_cor + fld_ind
      } else {
        fld = cbind(fld_cor, fld_ind)
      }
    }
  }
  
  attr(fld, "coefficients") = b
  attr(fld, "loadings") = Lambda
  attr(fld, "response.names") = object$response.id
  
  return(fld)
}

sp_fld = get_field_fixed(
  mB2,
  domain,
  which.response = species_idx
)

beta_fixed = fixef(mB2)$cond

species_levels = unique(as.character(dat$spid))
# Repeats the species ID to match the number of rows in 'domain', 
# then converts it to a categorical factor locked to the master 'species_levels' list.
# (Note: 'levels' does not change the repeated data
domain_j = data.frame(spid = factor(rep(target_species_id, nrow(domain)),
                                    levels = species_levels),
                      mi = domain$mi,
                      tempann = domain$tempann)

X_domain_j = model.matrix( ~0 +spid + spid:mi + spid:tempann, data = domain_j)

fixed = X_domain_j %*% beta_fixed
spatial = as.numeric(sp_fld[, 1])

domain$log_intensity = fixed + spatial

library(stars)
domain_stars = st_as_stars(domain[, c("lon", "lat", "log_intensity")], dims = c("lon", "lat"))
grid_pred_lonlat = domain_stars

species_points = dat[dat$spid == target_species_id & dat$occ == 1, ]
species_points_lonlat = st_as_sf(species_points, coords = c("lon", "lat"), crs = 4326)

border_lonlat = disBorder("NSW")

tmap_mode("plot")

load("m_intensity_map.RDATA")
mlgcp_map_multi_u = m_intensity_map(grid_pred_lonlat, model_label = "mlgcp, 2 l.v.")
print(mlgcp_map_multi_u)
save(mlgcp_map_multi_u, file = "results/mlgcp_map_multi_u.RDATA")


u = get_bfs_coeff(mB2)
lambda = get_loadings(mB2)

load("identify_ulambda.RDATA")
ULambda = identify_ulambda(u, lambda)
U_mlgcp = ULambda$U
Lambda_mlgcp = ULambda$Lambda
save(U_mlgcp, Lambda_mlgcp, file = "results/U_Lambda_mlgcp.RDATA")
