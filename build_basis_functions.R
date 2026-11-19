bfs <- make_basis(k = 50, domain, from.package = "scampr")

# Pruning basis functions with zero support
if (is(bfs, "scampr.bf")) {
  tmp <- bf_matrix(bfs, rbind(pres[,c("xm", "ym")], quad[,c("xm", "ym")])) # via the data?
  # tmp <- bf_matrix(bfs, pres[,c("xm", "ym")]) # via the data?
  bf.check <- apply(tmp, 2, function(x){sum(x > 0)})
  bfs <- bfs[bf.check != 0, ]
  rm(tmp, bf.check)
}

# Evaluate basis functions at data locations
Z_star = bf_matrix(bfs, dat[, c("x", "y")])
#Z_star = Z_star/max(apply(Z_star,2,sd)) #some standardisation needed for ridge fits
save(Z_star, file = "Z_star.RDATA")
load("Z_star.RDATA")