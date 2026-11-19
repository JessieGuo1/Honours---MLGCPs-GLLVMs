# Get the data
source("get_and_format_flora_data.R")
rm(dat_pa)

# Get a vector of the factor terms
facts <- c("disturb", "soilfert", "vegsis")

# Load the species information
load("flora_species.RDATA")

# Retain selected flora groups
flora_groups <- unique(sp.info$group)
pres <- pres[pres$group %in% flora_groups, ]

# Keep species with at least 10 presences
pres <- pres[pres$spid %in% names(table(pres$spid))[table(pres$spid) >= 10], ]

# Retained species
species <- unique(pres$spid)

# Dimensions
m <- length(species) # responses
q <- nrow(quad) # quadrature points

# Replicate quadrature points for each species
tmp <- quad[rep(seq_len(q), m), ]
tmp$spid <- rep(species, each = q) # assign the species id's

dat <- rbind(pres, tmp)
n <- nrow(dat)
rm(tmp)

# Preserve lon/lat and use model coordinates as x/y
dat$lon <- dat$x; dat$lat <- dat$y; dat$x <- dat$xm; dat$y <- dat$ym; dat$xm <- NULL; dat$ym <- NULL
save(dat, file = "dat.RDATA")
load("dat.RDATA")

load("domain.RDATA")
domain$lon <- domain$x; domain$lat <- domain$y; domain$x <- domain$xm; domain$y <- domain$ym; domain$xm <- NULL; domain$ym <- NULL

gc()