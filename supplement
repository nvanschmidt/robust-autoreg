# Code for autoregressive model. Includes conditional psi estimation,
# BRM calculation, IFM calculation, and BRM/IFM standardization across years. # This code was designed to work in conjunction with PRESENCE v11.0.

#########################################################################
############################### FUNCTIONS ###############################
#########################################################################

##### Import PRESENCE output #####

# Imports multiple "ocep" (reformatted PRESENCE model output) files into R.
# Imported ocep files must be named "ocep_i[iteration]_[metric][radius].csv"
# e.g., "ocep_i2_br5.csv") and be formatted in four columns:
# parameter (string psi, gam, eps, or P), season (numeric; integer except for
# P, which is season.visit), site (string), and estimate (parameter value).
# Note: this code uses decimals to order visits within a season and thus if
# number of of visits is >9 the columns will need to be manually reordered.
# Outputs a list of each radius, each of which contains a list of matrices of parameter values (o=psi,c=gamma,e=epsilon,p=p) 
import.ocep <- function (iteration, metric, r.list) {  # r.list is a numeric vector of radiuses used
  require(ecodist)
  ocep.list <- lapply(r.list, function(r) {
    ocep <- read.csv(paste("ocep_i",iteration,"_",metric,r,".csv",sep=""))
    ocep.split <- split(ocep, ocep$parameter)
    o <- crosstab(rowlab=site, collab=season, values=estimate,                                            data=ocep.split$psi)
    c <- cbind(NA,crosstab(rowlab=site, collab=season, values=estimate, data=ocep.split$gam))
    colnames(c) <- 1:ncol(c)
    e <- cbind(NA,crosstab(rowlab=site, collab=season, values=estimate, data=ocep.split$eps))
    colnames(e) <- 1:ncol(e)
    p <- crosstab(rowlab=site, collab=season, values=estimate, data=ocep.split$P)
    list(o=o,c=c,e=e,p=p)
  })
  names(ocep.list) <- paste(paste(metric,r.list,sep=""))
  return(ocep.list)
}

# Performs an alignment check against the survey history, ocep list, and any
# other arguments to ensure they imported correctly and that site names match.
check.alignment.list <- function (dh, ocep.list, ...) {
  lapply(ocep.list, function (x) {
    sapply(list(x$o, x$c, x$e, x$p, ...), function (x) { all.equal(rownames(dh),rownames(x)) })
    })}


##### Calculate Conditional Psi ##### 

# Calculates a "survey history" (1 = detected in that season, 0 = not detected in that season or
# unsurveyed) from a detection history (1 = detection, 0 = non-detection, NA = not surveyed).
# Used in function get.az.
get.sh <- function(dh,num.visits,season.list) {
  dh.0 <- dh
  dh.0[is.na(dh.0)] <- 0
  output <- as.data.frame(t(sapply(1:nrow(dh.0), function(i) {
    sapply(1:length(num.visits), function(k) {
      max(dh.0[i,(1+sum(num.visits[0:(k-1)])) : sum(num.visits[1:k]) ])
    })})),row.names=rownames(dh))
  colnames(output) <- season.list
  return(output)
}

# Gets possible alternate Z (latent unknown occupancy histories).
# Excludes Z with 0 probability (species not present when detected).
get.az <- function(sh) {
  all.az <- data.matrix(expand.grid(rep(list(0:1),ncol(sh))))
  apply(sh, MARGIN = 1, FUN = function(x) unique(t(pmax(t(all.az),x))))
}

# Parallel processing calculates conditional psi based on alternate Z possibilites.
get.cpsi <- function (az, o, c, e, p, nv, unsurveyed=NULL) { 
  # creates local copies of inputs with only surveyed sites
  # (saves computation time, o = cpsi for unsurveyed sites)
  o.s <- o[!(row.names(o) %in% unsurveyed), ]
  o.u <- o[row.names(o) %in% unsurveyed, ]
  c.s <- c[!(row.names(c) %in% unsurveyed), ]
  e.s <- e[!(row.names(e) %in% unsurveyed), ]
  p.s <- p[!(row.names(p) %in% unsurveyed), ]
# creates local copies of remaining data for use in parallel processing
  az.s <- az
  nv.s <- nv
  cl <- makeCluster(nc) # make parallel processing cluster based on number CPU cores
  azp.all <- parLapply(cl=cl, 1:length(az.s), function(i) {  # loop over sites
    # create empty product (transition * detection, per season) probability matrix
    pp <- matrix(NA,nrow=nrow(az.s[[i]]),ncol=ncol(az.s[[i]]))
    # loop over rows of alternate possible true occupancy states [z] to to fill in probabilities
    sapply(1:nrow(az.s[[i]]), function(z) {  
      ifelse(az.s[[i]][z,1] == 1,   # fill in first column for psi; if occupied:
            {
             # Probability of not being detected over all visits in that season.
             # Also used for years when species detected for coding convenience
             # because it cancels out in those cases regardless.
             p.nd <- 1; for (k in 1:nv.s[1]) { p.nd <- p.nd*(1-p.s[i,k]) }
             pp[z,1] <<- o.s[i,1]*p.nd # extinction probability * detection probability
            },
             pp[z,1] <<- 1-o.s[i,1])  # else probability unoccupied
      sapply(2:ncol(az.s[[i]]), function(j) {  # loop over rest of seasons [j] to fill in epsilon/gamma
        ifelse(az.s[[i]][z,j-1] == 1,  # if occupied in previous season:
               ifelse(az.s[[i]][z,j] == 1,  # if occupied in current season:
                      {
                      p.nd <- 1; past.k <- sum(nv.s[1:j-1]); for (k in 1:nv.s[j]) {
                        p.nd <- p.nd*(1-p.s[i,past.k+k]) } #non-detection probability; = 1 if not surveyed
                        pp[z,j] <<- (1-e.s[i,j])*p.nd # extinction probability * non-detection probability
                      },
                      pp[z,j] <<- e.s[i,j]),  # else probability went extinct
                                     # else unoccupied in previous season
               ifelse(az.s[[i]][z,j] == 1,  # if occupied in current season:
                      {               
                      p.nd <- 1; past.k <- sum(nv.s[1:j-1]); for (k in 1:nv.s[j]) {
                        p.nd <- p.nd*(1-p.s[i,past.k+k]) }
                        pp[z,j] <<- c.s[i,j]*p.nd  # colonization probability * detection probability 
                      },
                      pp[z,j] <<- 1-c.s[i,j]))  # else probability not colonized
      })})
    ppp <- apply(pp, MARGIN = 1, FUN = prod)  # overall probability of each possible Z
    # probabilities of alternate Zs for this site where it was occupied 
    # in that season; appends ppp for calculating cpsi below
    azp <- cbind(az.s[[i]]*ppp, ppp)
    return(azp)})
  stopCluster(cl)  # stops parallel processing
  cpsi <- matrix(unlist(lapply(azp.all, function(i) {  # calculates conditional psi
    sum.ppp <- sum(i[,'ppp'])
    # sum probability of Z where site was occupied in that season,
    # divided by sum probability of all Z for site (Eq. 4 in text)
    apply(i[,1:ncol(i)-1,drop=FALSE], MARGIN = 2, FUN = function(j) sum(j)/sum.ppp)
  })),byrow=TRUE,nrow=length(azp.all))
  row.names(cpsi) <- row.names(o.s)  # retrieves names of surveyed sites for output table
  cpsi <- rbind(cpsi, o.u)  # appends unsurveyed sites' unconditional occupancy 
  cpsi <- cpsi[order(rownames(cpsi)),]  # reorders output to match input (alphabetizes)
  colnames(cpsi) <- paste("CPSI.",season.list,sep="")
  return(cpsi)
}

# Runs the above function for a list of different ocep tables
# (e.g., from different radii or other alternative parameterizations). 
get.cpsi.list <- function(az, ocep.list, nv, unsurveyed=NULL) {
  cpsi.list <- lapply(ocep.list, function (x) {
    get.cpsi(az,x$o,x$c,x$e,x$p,nv, unsurveyed)
  })
  names(cpsi.list) <- names(ocep.list)
  return(cpsi.list)
}


##### Connectivity Metric Calculations #####

# Calculates buffer radius metrics using conditional psi.
# Takes as input conditional psi matrix, pairwise distance matrix
# between sites, matrix of patch areas in each season, and radius.
get.brm <- function(cpsi,dist,area,r) {
  brm <- cpsi  # creates matrix to fill in loop
  dist[dist <= r] <- 1   # includes sites within radius
  dist[dist > r] <- 0  # excludes sites outside of radius
  dist <- dist * (1-diag(nrow(dist)))  # excludes focal site
  for (j in 1:ncol(cpsi)) {  # loops over seasons
    for (i in 1:nrow(cpsi))  {  # loops over sites
      brm[i,j] <- sum(cpsi[,j] * area[,j] * dist[i,])  # calculates metric
    }}
  return(brm)
}


# Calculates buffer radius metrics (log10 + 1) for multiple radii and
# corresponding conditional psi values from that radius' autoregressive model
# (e.g., radius = 1 for cpsi estimated from a radius 1 autoregressive model;
# radius 2 for cpsi estimated from a radius 2 autoregressive model, etc.).
# Requires cpsi naming convention: cpsi.i[iteration].br$br[radius]
get.brm.list <- function (r.list,cpsi.list,season.list,iteration) {
  brm.list <- lapply(r.list, function(r) {
    brm <- get.brm(eval(parse(text=gsub("T",r,paste("cpsi.i",iteration,".br$brT",sep="")))), vas.dist, vas.area, r)
    lbrm <- log10(brm + 1)
    return(lbrm)
  })
  output <- matrix(unlist(brm.list), ncol = length(r.list) * ncol(cpsi.list[[1]]), byrow = FALSE)  # reorganizes output into matrix
  rownames(output) <- row.names(cpsi.list[[1]])  # assigns row headings (patch name)
  colnames(output) <- unlist(lapply(r.list, function (r) {paste("LB",r,".",season.list,sep="")}))  # assigns column headings (LB[radius].[season])
  return(output)
}

# As above function (but for the first iteration when there is only one model estimating conditional psi)
get.brm.list.i1 <- function (r.list,cpsi,season.list) {
  brm.list <- lapply(r.list, function(r) {
    brm <- get.brm(cpsi, vas.dist, vas.area, r)
    lbrm <- log10(brm + 1)
    return(lbrm)
  })
  output <- matrix(unlist(brm.list), ncol = length(r.list) * ncol(cpsi), byrow = FALSE)
  rownames(output) <- row.names(cpsi)
  colnames(output) <- unlist(lapply(r.list, function (r) {paste("LB",r,".",season.list,sep="")}))
  return(output)
  }

# Calculates incidence function metrics using conditional psi.
# Takes as input conditional psi matrix, pairwise distance matrix
# between sites, matrix of patch areas in each season, and "radius" (1/alpha).
get.ifm <- function(cpsi,dist,area,r)  {
  ifm <- cpsi  # creates matrix to fill in loop
  not.self <- 1-diag(nrow(dist))  # excludes focal site from calculation
  for (j in 1:ncol(cpsi)) {  # loops over seasons
    for (i in 1:nrow(cpsi))  {  # loops over sites
      ifm[i,j] <- sum(not.self[i,] * cpsi[,j] * area[,j] * exp(-1/r * dist[i,]))  #calculates metric
    }}
  return(ifm);
}

# Calculates incidence function metric (log10 + 1) for multiple "radii" and
# corresponding conditional psi values from that "radius"' autoregressive model.
# Requires cpsi naming convention: cpsi.i[iteration].if$if[radius]
get.ifm.list <- function (r.list,cpsi.list,season.list,iteration) {
  ifm.list <- lapply(r.list, function(r) {
    ifm <- get.ifm(eval(parse(text=gsub("T",r,paste("cpsi.i",iteration,".if$ifT",sep="")))), vas.dist, vas.area, r)
    lifm <- log10(ifm + 1)
    return(lifm)
  })
  # reorganizes output into matrix
  output <- matrix(unlist(ifm.list), ncol = length(r.list) * ncol(cpsi.list[[1]]), byrow = FALSE)
  rownames(output) <- row.names(cpsi.list[[1]])  # assigns row headings (patch name)
  # assigns column headings (IF[radius].[season])
  colnames(output) <- unlist(lapply(r.list, function (r) {paste("IF",r,".",season.list,sep="")}))  
  return(output)
}

# As above function, but for the first iteration (when there is only one model estimating conditional psi)
get.ifm.list.i1 <- function (r.list,cpsi,season.list) {
  ifm.list <- lapply(r.list, function(r) {
    ifm <- get.ifm(cpsi, vas.dist, vas.area, r)
    lifm <- log10(ifm + 1)
    return(lifm)
  })
  output <- matrix(unlist(ifm.list), ncol = length(r.list) * ncol(cpsi), byrow = FALSE)
  rownames(output) <- row.names(cpsi)
  colnames(output) <- unlist(lapply(r.list, function (r) {paste("IF",r,".",season.list,sep="")}))
  return(output)
}

# Standardizes buffer radius or incidence function metrics for a
# given radius, across seasons but excluding the final season
# because it is not used in modelling (metrics are lagged one season as gamma & epsilon covariates).
standardize.metric <- function(metric,r.list,season.list) {
  rl <- length(r.list)
  yl <- length(season.list)
  means <- sapply(0:(rl-1), function (r) {  # loops over radii
    mean <- mean(metric[,(1+r*yl):(yl-1+r*yl)])  # calculates mean across seasons (except final season)
    rep(mean,yl)  # fills in mean across seasons
  })
  sds <- sapply(0:(rl-1), function (r) {  # loops over radii
    sd <- sd(metric[,(1+r*yl):(yl-1+r*yl)]) # calculates st. dev. across seasons (except final season)
    rep(sd,yl)  # fills in standard deviation across seasons
  })
  scale(metric,center=c(means),scale=c(sds))  # standardizes metrics
}


########################################################################
#################### EXAMPLE OF RUNNING THE MODEL ######################
########################################################################

##### Setup: load parallel processing #####
require(parallel)
nc<-detectCores()

##### Setup: load core data inputs #####

# Define number of visits in each season, season names, and radius (or 1/alpha) values
num.visits <- c(5,rep(3,13)) # number of visits in each season
season.list <- c("02", "03", "04","05","06","07","08","09","10","11","12","13","14","15")
r.list <- c(1:10, 15, 20, 25, 30)

# Load & alphabetize detection history
# Expected format: rows = patches, columns = visits,
# first row = visit names, first column = patch names (SITE)
# Expected values: 1 = detection, 0 = non-detection, NA = not surveyed
dh <- read.csv("blra_dh.csv",row.names="SITE")
dh <- dh[order(rownames(dh)),]

# Load & alphabetize pairwise distance matrix (diagonal identity = 0)
dist <- as.matrix(read.csv("patch_distances.csv",row.names="SITE"))
dist <- dist[order(rownames(dist)),]
dist <- dist/1000 # rescale from meters to kilometers

# Load & alphabetize patch areas (annual measure; one column per season)
area <- read.csv("patch_areas.csv", row.names="SITE")
area <- area[order(rownames(area)),]
area <- area/10000 # rescale square meters to hectares

# Get list of unsurveyed sites
unsurveyed <- row.names(dh[rowSums(is.na(dh)) == ncol(dh),])

# Get survey history and alternate Z possibilities for surveyed sites
sh <- get.sh(dh[!(row.names(dh) %in% unsurveyed), ],num.visits,season.list)
az <- get.az(sh)


########## Iteration 1 ###########

### Prior to CPSI/BRM/IFM calculation the model must be fit in Program PRESENCE using
### a multi-season single species occupancy model, with all unsurveyed sites included
### (with detection histories of missing data for all visits and all appropriate site
### covariate information). The parameter estimates in the model's .out file must then
### be reformatted into an OCEP table (see comments for function import.ocep() above).

# Import OCEP file (must be in current working directory) and check alignment.
ocep.i1 <- import.ocep(1,"br",0)
check.alignment.list(dh=dh,ocep.list=ocep.i1,area,dist)
# Make sure all values from above checks are TRUE before proceeding.
# If any are FALSE, re-check formatting and import of OCEP files.

# Calculate conditional psi
cpsi.i1.br <- get.cpsi.list(az,ocep.i1,num.visits,unsurveyed)

# Calculate, standardize, and export BRMs
brm.i1 <- get.brm.list.i1(r.list,cpsi.i1.br$br0,season.list)
brm.i1.s <- standardize.metric(brm.i1,r.list,season.list)
write.csv(brm.i1.s,"brm_i1_s.csv")

# Calculate, standardize, and export IFMs
ifm.i1 <- get.ifm.list.i1(r.list,cpsi.i1.br$br0,season.list)
ifm.i1.s <- standardize.metric(ifm.i1,r.list,season.list)
write.csv(ifm.i1.s,"ifm_i1_s.csv")

########## Iteration 2 ########### 

### The BRM/IFMs calculated above must be entered into a new PRESENCE project
### and a new modelset must be run, adding them as covariates on psi, gamma,
### and/or epsilon. The parameter estimates in the .out file should then be
### reformatted into an ocep table, as before. If running multiple radius values
### each model (i.e., each radius) should be fit to as a separate model and
### imported as a separate ocep table file.

# Import OCEP file (must be in current working directory) and check alignment.
ocep.i2.br <- import.ocep(2,"br",r.list)
check.alignment.list(dh=dh,ocep.list=ocep.i2.br,area,dist)
ocep.i2.if <- import.ocep(2,"if",r.list)
check.alignment.list(sh=sh,ocep.list=ocep.i2.if,area,dist)
# Make sure all values from above checks are TRUE before proceeding.
# If any are FALSE, re-check formatting and import of OCEP files.

# Calculate conditional psi
cpsi.i2.br <- get.cpsi.list(az,ocep.i2.br,unsurveyed)
cpsi.i2.if <- get.cpsi.list(az,ocep.i2.if,unsurveyed)

# Calculate, standardize, and export BRMs
brm.i2 <- get.brm.list(r.list,cpsi.i2.br,season.list,2)
brm.i2.s <- standardize.metric(brm.i2,r.list,season.list)
write.csv(brm.i2.s,"brm_i2_s.csv")

# Calculate, standardize, and export IFMs
ifm.i2 <- get.ifm.list(r.list,cpsi.i2.if,season.list,2)
ifm.i2.s <- standardize.metric(ifm.i2,r.list,season.list)
write.csv(ifm.i2.s,"ifm_i2_s.csv")

### Repeat iteration process (copy and paste Iteration 2 code, changing
### the iteration numbers each time) until parameter values of all radius'
### models converge (generally ~10 iterations). 
