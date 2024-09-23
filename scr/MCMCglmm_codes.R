# R codes used in the study to fit the MCMCglmms models are available below.

library(MCMCglmm)
# This script shows as a sample the code to fit the MCMCglmms models.

# MCMCglmm analyses
## Binomial models
priors <- list(R = list(V = 1, nu = 50),
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

# Read the phylogenetic tree ####
phangorn::nnls.tree(cophenetic(ape::read.tree("data/seedArc_treeV0.tree")),
                    ape::read.tree("data/seedArc_treeV0.tree"), method = "ultrametric") -> nnls
nnls$node.label <- NULL
nnls.a<- drop.tip(nnls, unique(data$animal))
nnls<- drop.tip(nnls, nnls.a$tip.label)

nite = 500000
nthi = 100
nbur = 100000

# Models ####
model.par <-c(
  "scale(Tmean)+scale(Tmean):scale(p1)","scale(Tmean)+scale(Tmean):scale(p2b)",
  "scale(Alternating)+scale(Alternating):scale(p1)","scale(Alternating)+scale(Alternating):scale(p2b)",
  "scale(Scarification)+scale(Scarification):scale(p1)","scale(Scarification)+scale(Scarification):scale(p2b)",
  "scale(cold_str)+scale(cold_str):scale(p1)","scale(cold_str)+scale(cold_str):scale(p2b)",
  "scale(warm_str)+scale(warm_str):scale(p1)","scale(warm_str)+scale(warm_str):scale(p2b)",
  "scale(Fire)+scale(Fire):scale(p1)","scale(Fire)+scale(Fire):scale(p2b)")

### List of random factors ####
random_factors <- list(
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"),
  c("animal", "ID", "id_test", "doi", "seedlot", "substrate"))

# RUN MODELS global p1 and p2  ####
library(doParallel); library(foreach)

parallel::detectCores()
n.cores <- 2
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

clusterEvalQ(my.cluster, {
  library(MCMCglmm)
})

clusterExport(my.cluster, varlist = c("data", "nnls", "nite", "nthi", "nbur", "priors", "model.par", "random_factors"))
model.res <- foreach(i = 1:length(model.par), .packages = "foreach") %dopar% {
  fixed <- as.formula(paste("cbind(Germinated, Germinable-Germinated) ~", model.par[i], sep = ""))
  random <- random_factors[[i]]
  
  mm <- MCMCglmm(fixed = fixed, 
                 random = as.formula(paste("~", paste(random, collapse = "+"))),
                 family = "multinomial2", pedigree = nnls, prior = priors, data = data,
                 nitt = nite,thin = nthi,burnin = nbur,verbose = FALSE
  )
}
parallel::stopCluster(cl = my.cluster)
