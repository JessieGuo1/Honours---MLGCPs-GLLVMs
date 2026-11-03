# set approximate number of basis function nodes to use
basis.dim = 50

library(glmmTMB)

# get the data
source("get_and_format_flora_data.R")
rm(dat_pa)

# get the relevant predictors according to pred_set (these are all basically derived from elevation, we are just looking at moisture and temp. covars and soil fertility score and disturbance index)
# preds <- switch(pred_set,
#                 A = "None",
#                 B = c("mi", "tempann"),
#                 C = c("mi", "tempann", "mi2", "tempann2"),
# )

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
# Fit an intercepts-only model #################################################
################################################################################

# need a function to strip environments from the model objects so they can be saved without being massive
#rm.env <- function(m) {
#  attr(attr(m$frame, "terms"), ".Environment") <- NULL
#  attr(m$modelInfo$terms$cond$fixed, ".Environment") <- NULL
#  attr(m$modelInfo$reTrms$cond$terms$fixed, ".Environment") <- NULL
#  attr(m$modelInfo$allForm$combForm, ".Environment") <- NULL
#  attr(m$modelInfo$allForm$formula, ".Environment") <- NULL
#  attr(m$modelInfo$allForm$ziformula, ".Environment") <- NULL
#  attr(m$modelInfo$allForm$dispformula, ".Environment") <- NULL
#  m
#}

################################################################################
# Fit the model including linear env terms #####################################
################################################################################

# with correlated fields only ##################################################
time_i = system.time({try(assign("mB1", mvlgcp(occ ~ 0 + spid + spid:mi + spid:tempann,
                         data = dat, weights = wt, basis.functions = bfs, response.id = spid, n_factors = 1,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))})

# save if the model fitted
if (exists("mB1")) {
  print("mB1 converged without start parameters (attempt 1)")
#  mB1 <- rm.env(mB1)
  mB1$start.pars <- "cold"
  save(mB1, file = paste0(getwd(), "/models/mB1_1.RDATA"))
}

load("models/mB1_1.RDATA")
mB1

pred_mlgcp = predict(mB1, type = "response", re.form = NULL)

loglik_new_mlgcp = sum(dat$occ * log(pred_mlgcp) - dat$wt * pred_mlgcp)
save(loglik_new_mlgcp, file = "loglik_new_mlgcp.RDATA")

time_mlgcp = as.numeric(time_i[3])
time_mlgcp
save(time_mlgcp, file = "time_mlgcp.RDATA")

##########################
#Using tmap
target_species_id <- "nsw19"
Lambda = get_loadings(mB1)
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
  mB1,
  domain,
  which.response = species_idx
)

beta_fixed = fixef(mB1)$cond

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
mlgcp_map = m_intensity_map(grid_pred_lonlat, species_label = 19)
print(mlgcp_map)
save(mlgcp_map, file = "mlgcp_map.RDATA")
load("mlgcp_map.RDATA")
print(mlgcp_map)










#Inhomogeneous K func
# Diagnostic for chosen species, target_species_id
# Calculate observed, simulated inhomogeneous K functions
# Select species

domain_j$response.id = domain_j$spid
# Propagate variance of spatially varying term through smoothed terms at every grid location
eta_full = predict(mB1, newdata = domain_j, type = "link", re.form = NULL) #re.form controls whether random effects included, NULL = include
eta_fixed = predict(mB1, newdata = domain_j, type = "link", re.form = NA) 
spatial = eta_full - eta_fixed
var_spatial = var(spatial)

# This species' fixed-effect linear predictor - used as mean
domain$mean = fixed

domain$intensity = exp(domain$log_intensity)

# Build observation window in same coordinate system used throughout
border_xy = st_transform(border_lonlat, crs = 28356)
# Convert projected coordinates to kilometres to match domain$x/y/xm/ym
border_geom = st_union(st_geometry(border_xy)) / 1000
# Remove CRS metadata for spatstat
st_crs(border_geom) = NA
win = as.owin(st_sf(geometry = border_geom))

# Convert fixed-effect log-intensity surface and full fitted intensity surface into spatstat pixel images
fixed.im = as.im(domain[, c("x", "y", "mean")], W = win)
pred.im = as.im(domain[, c("x","y","intensity")], W = win)

# Rename pres_df's spid to match dat_rename's convention, then subset to just this species' own presence locations
pres_df_rename = pres %>% 
  mutate(spid = as.numeric(str_extract(spid, "\\d+")))
sp_presence_df = pres_df_rename[pres_df_rename$spid == target_species_id, ]

# Build species' point pattern within study window using same kilometre-scale as win
po_sp_ppp = ppp(x = sp_presence_df$xm,
                y = sp_presence_df$ym,
                window = win)

#load("opt.RDATA")

# Observed inhomogeneous K-function for this species, using its fitted intensity surface
mlgcp_K_obs <- Kinhom(po_sp_ppp, lambda = pred.im, correction = "border")

# Simulate comparison point patterns from a Matern log-Gaussian Cox process and construct simulation envelope
mlgcp_K_env <- envelope(po_sp_ppp, fun = Kinhom, simulate = expression(rLGCP(model = "matern", mu = fixed.im, var = var_smooth, nu = 3/2, scale = bfs$scale[1])))

# Plot observed K-function and simulation envelope
load("Kfunc.RDATA")
Kfunc(mlgcp_K_obs, mlgcp_K_env)
mlgcp_kfunc = recordPlot()
save(mlgcp_kfunc, file = "mlgcp_kfunc.RDATA")

load("mlgcp_kfunc.RDATA")
replayPlot(mlgcp_kfunc)




# Fitted Fields ################################################################
home.wd <- getwd()

selected.sp <- c(2, 4, 5, 6, 7, 14)
source("get_and_format_flora_data.R")
domain$lon <- domain$x; domain$lat <- domain$y; domain$x <- domain$xm; domain$y <- domain$ym; domain$xm <- NULL; domain$ym <- NULL
rm(dat_pa)
load("flora_species.RDATA")

# get the groups to which the flora species belongs to subset the data
flora_groups <- unique(sp.info$group)

# subset the presence-only data (in the case we are using restricted species)
pres <- pres[pres$group %in% flora_groups, ]

# want to use species with at least 10 presences
pres <- pres[pres$spid %in% names(table(pres$spid))[table(pres$spid) >= 10], ]

species <- unique(pres$spid)

# re-order the species information
sp.info <- sp.info[order(sp.info$spid), ]
tmp <- strsplit(sp.info$full_name, " ", fixed = T)
sp.info$nice.names <- paste(substr(sapply(tmp, function(x){x[1]}), start = 1, stop = 1), substr(sapply(tmp, function(x){x[2]}), start = 1, stop = 4), sep = ". ")
rm(tmp)
# need to fix the duplicate names of one species (likely an error from the paper)
sp.info$nice.names[sp.info$nice.names == "S. lueh"] <- paste0(sp.info$nice.names[sp.info$nice.names == "S. lueh"], c(" A", " B"))

# crop the sp.info
sp.info <- sp.info[sp.info$spid %in% unique(pres$spid), ]

sp_flds <- get_field(mB1, domain, which.response = selected.sp)
tmp.flds <- get_field(mB1, domain)
lat_flds <- tmp.flds[ , c("d1", "d2")]

zlims <- range(cbind(sp_flds, lat_flds))
tmp.lamb <- get_loadings(mB1)
Lambda <- tmp.lamb[selected.sp, ]
rm(tmp.flds, tmp.lamb)

#tmp <- VarCorr(mB1)
#sigma_j <- diag(tmp$cond$basis.functions.1)[selected.sp]
#sigma_j <- c("3.2e-8", "2.6e-9", "2.6e-10", "1.1e-9", "2.1e-8", "2.007")
dir.create(file.path(home.wd, "figures"))

plot.res <- 500
png(filename = paste0(home.wd, "/figures/nsw_fitted_fields.png"), width = 6.2 * plot.res, height = 3.5 * plot.res, res = plot.res)
layout(matrix(c(1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,13,14,15,16,17,18), nrow = 3, ncol = 7, byrow = T), widths = c(0.1, rep(0.9/6, 6)), heights = rep(1/3,3))
par(mar = c(0,0.1,1.5,0))

plot(1, axes = F, bty = "n", type = "n")

plot(vec2im(sp_flds[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(pres[pres$spid == sp.info$spid[selected.sp][1], c("xm", "ym")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = sp.info$nice.names[selected.sp][1])
mtext(expression(paste(ln~mu[j],"(s)")), side = 2, line = 0.5, las = 1, xpd = T)
mtext(bquote(lambda[1]==.(round(Lambda[1,1],1))~","~.(round(Lambda[1,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(pres[pres$spid == sp.info$spid[selected.sp][2], c("xm", "ym")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = sp.info$nice.names[selected.sp][2])
mtext(bquote(lambda[2]==.(round(Lambda[2,1],1))~","~.(round(Lambda[2,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)


plot(vec2im(sp_flds[,3], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(pres[pres$spid == sp.info$spid[selected.sp][3], c("xm", "ym")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = sp.info$nice.names[selected.sp][3])
mtext(bquote(lambda[3]==.(round(Lambda[3,1],1))~","~.(round(Lambda[3,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,4], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(pres[pres$spid == sp.info$spid[selected.sp][4], c("xm", "ym")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = sp.info$nice.names[selected.sp][4])
mtext(bquote(lambda[4]==.(round(Lambda[4,1],1))~","~.(round(Lambda[4,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,5], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(pres[pres$spid == sp.info$spid[selected.sp][5], c("xm", "ym")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = sp.info$nice.names[selected.sp][5])
mtext(bquote(lambda[5]==.(round(Lambda[5,1],1))~","~.(round(Lambda[5,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,6], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(pres[pres$spid == sp.info$spid[selected.sp][6], c("xm", "ym")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = sp.info$nice.names[selected.sp][6])
mtext(bquote(lambda[6]==.(round(Lambda[6,1],1))~","~.(round(Lambda[6,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(1, axes = F, bty = "n", type = "n")

plot(vec2im(lat_flds[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(xi[1]^(u),"(s)")), side = 2, line = -3, las = 1, xpd = T)

plot(vec2im(lat_flds[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(xi[2]^(u),"(s)")), side = 2, line = -3, las = 1, xpd = T)

par(mar = c(0,2,4,1))
legend_image <- as.raster(matrix(rev(viridisLite::plasma(20)), ncol=1))
plot(c(0,8),c(zlims[1],zlims[2]),type = 'n', axes = F, xlab = "", ylab = "", main = "")
mtext("field values", line = 0.5, cex = 0.75, adj = 0, at = 1.75)
text(x=1.9 + 3, y = seq(round(zlims[1]),round(zlims[2]),l=5), labels = seq(round(zlims[1]),round(zlims[2]),l=5), xpd = T)
lines(c(0,1.01,1.01,0,0) + 3, c(zlims[1],zlims[1],zlims[2],zlims[2],zlims[1]))
for (i in seq(round(zlims[1]),round(zlims[2]),l=5)) {
  lines(c(1.01, 1.21) + 3, rep(i, 2))
}
rasterImage(legend_image, 0 + 3, zlims[1], 1 + 3, zlims[2])

par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

