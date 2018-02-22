###################################################################
# Species-Habitat Realationships of Urban Bats at Multiple Scales #
################ By Travis Gallo and Mason Fidino #################
#################### Urban Wildlife Institute #####################

# Function to load and download packages (if needed)
package_load <- function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # Download required packages if they're not already
  
  pkgsToDownload <- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # Dhen load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# The packages needed for this analysis
packs <- c("dplyr", "reshape2", "runjags", "mcmcplots",
          'parallel', 'raster','rgdal','stringr',
          'foreach','doParallel')

package_load(packs)

# Read in bat detection data
sp_dat <- readRDS("gallo_et_al_2018_bat_detections.RDS")

# We will be doing a lot of nested indexing within the model
# This function makes two indexes
make_index <- function(x){
  yrses <- as.numeric(factor(x$year.session)) # survey year
  sites <- as.numeric(factor(x$site)) # sampling site
  return(data.frame(yrs = yrses, sites = sites))
}

# We use the data for EPFU here, but this index will work for all three species
indx <- make_index(sp_dat$epfu)

# Make a final index that indicates which year a session was sampled
yrmu_vec <- as.numeric(factor(substr(unique(sp_dat$epfu$year.session), 1,4)))

# Read in how many days each site was sampled during each session (j-matrix)
jmatlong <- readRDS("gallo_et_al_2018_j_matrix.RDS")

# Determine sites where no observations occured
to_na <- which(jmatlong$value==0)

# Detection data needs to have an NA value when no observations occured
# Add associated NA values to sp_dat
for(i in 1:length(sp_dat)){
  sp_dat[[i]]$spcount[to_na] <- NA
}

# load landcover covariates at the three scales
# R script that extracts these landcovers is avaliable at https://github.com/TravisGallo/Functions
# For the sake of privacy we are not disclousing coordinates for each site

orig_covs <- readRDS("gallo_et_al_2018_landcover_covariates.RDS")

# Scale the continuous covariates
aa <- orig_covs[,c(2:10)]
aa <- aa[,order(colnames(aa))]
aa <- scale(aa)
covs <- data.frame(site = orig_covs$site,aa,trt = as.character(orig_covs$trt), stringsAsFactors = FALSE)

# Make trt a binary variable
covs$trt[covs$trt=="U"] <- 1
covs$trt[covs$trt=="R"] <- 0
covs$trt <- as.numeric(covs$trt)

# Make site a factor
covs$site <- factor(covs$site)

# Create covariate matrix
x <- as.matrix(covs[,-1])

# Indicate the number of covariates we are using for JAGS model
ncov <- ncol(x)

# Function to run JAGs model for each species and extract parameters
model_run <- function(species) {
  # Choose the detection data for the species to be analyzed (epfu, labo, lano)
  choice <- sp_dat[[species]]$spcount

  # JAGS model data list
  data_list <- list(y = choice, x = x,
                  yrs = indx$yrs, sites = indx$sites, yr_mu_vec =yrmu_vec,
                  ncov = ncov, nyr = 3, nsession = length(yrmu_vec),
                  n = length(choice), jmat =  jmatlong$value)

  # Change the initial values for the z-parameter for each species
  z2 <- choice
  z2[z2>1] <- 1
  z2[is.na(z2)] <- 0
  
  # Create initial values with random starting values for EVERY parameter
  inits_fnc <- function(chain){
    gen_list <- function(chain = chain){
      list( 
        z = z2,
        laplace = runif(ncov, -3, 3),
        int_mu = runif(1, -3,3),
        int_yr_mu = runif(3, -3, 3),
        lp_mu = runif(1, -3, 3),
        lp_mu_yr = runif(3, -3, 3),
        int_tau_session = rgamma(1,1,1),
        int_tau_year = rgamma(1,1,1),
        lp_tau_session = rgamma(1,1,1),
        lp_tau_year = rgamma(1,1,1),
        lp_0 = runif(8, -3, 3),
        beta_0 = runif(8, -3, 3),
        pi = round(rbeta(ncov,1,1)),
        lambda = rgamma(1,1,1),
        pp=rbeta(1,1,1),
        
        .RNG.name = switch(chain,
                           "1" = "base::Wichmann-Hill",
                           "2" = "base::Marsaglia-Multicarry",
                           "3" = "base::Super-Duper",
                           "4" = "base::Mersenne-Twister",
                           "5" = "base::Wichmann-Hill",
                           "6" = "base::Marsaglia-Multicarry",
                           "7" = "base::Super-Duper",
                           "8" = "base::Mersenne-Twister"),
        .RNG.seed = sample(1:1e+06, 1)
      )
    }
    
    return(switch(chain,           
                  "1" = gen_list(chain),
                  "2" = gen_list(chain),
                  "3" = gen_list(chain),
                  "4" = gen_list(chain),
                  "5" = gen_list(chain),
                  "6" = gen_list(chain),
                  "7" = gen_list(chain),
                  "8" = gen_list(chain)
    )
    )
  }

  # JAGS model analysis
  model <- run.jags(model = "gallo_et_al_20198_bat_model_lasso.R",
                        monitor = c("beta", "lambda", "beta_0", "int_mu",
                                   "pi", "int_sd_session", "int_sd_year", "lp_0",
                                   "lp_mu", "lp_mu_yr"),
                        inits = inits_fnc,
                        data = data_list,
                        n.chains = 7, adapt = 5000, burnin = 20000, sample = 20000,
                        summarise = TRUE,
                        plots = FALSE,
                        method = "parallel")

  # Make the output a matrix
  md <- as.matrix(as.mcmc.list(model), chains = TRUE)
  
  # Estimate the probability each parameter should be included in model - "variable inclusion probability"
  best_pi <- md[, grep("pi\\[", colnames(md))] # extract Bernoulli results for each parameter
  par_prob <- apply(best_pi, 2, mean) # proportion of MCMC samples that got a 1
  names(par_prob) <- colnames(x) # parameter name as column name
  par_prob <- t(t(par_prob)) # transpose object
  
  # Selecting parameters to include in the model
  # Relative parameter selection, 
  # If the probability of the parameter being selected is < the mean probability of all parameters its centered on zero
  tozero <- which(par_prob<mean(par_prob)) # paramters that are less than mean
  mean.vip <- mean(par_prob) # mean vip
  
  # Get 95 % credible intervals
  best <- md[, grep("beta\\[", colnames(md))] # extract beta values
  # Set insignificant variables to 0
  for(i in 1:ncol(best)){
    zero <- which(best[,i]==0)
    if(i %in% tozero){
      best[-zero,i] <- best[-zero,i] - median(best[-zero,i])
    }
  }
  
  # Median and 95% confidence intervals
  # Create objects to store values
  o95 <- matrix(0, ncol = 3, nrow = c(ncol(best)))
  o90 <- matrix(0, ncol = 3, nrow = c(ncol(best)))
  m1 <- rep(0, ncol(best))
  
  # Loop through each beta parameter and estimate values
  for(i in 1:ncol(best)){
    n0 <- which(best[,i] != 0)
    m1[i] <- length(n0)
    o95[i,] <- quantile(best[n0,i], probs = c(0.025, 0.5, 0.975))
    o90[i,] <- quantile(best[n0,i], probs = c(0.05, 0.5, 0.95))
  }
  
  # Give matrix rownames
  rownames(o95) <- colnames(x)
  rownames(o90) <- colnames(x)
  colnames(best) <- colnames(x)
  
  # Reorder matrix for plotting
  o95 <- o95[c(10,7,9,8,4,6,5,1,3,2),]
  o90 <- o90[c(10,7,9,8,4,6,5,1,3,2),]
  par_prob <- par_prob[c(10,7,9,8,4,6,5,1,3,2)]
  best <- best[,c(10,7,9,8,4,6,5,1,3,2)]
  
  return(list(best=best, mean.vip=mean.vip, parprob=round(par_prob,2),o95=o95, m1=m1, o90=o90, model_out=model,model_mat=md))
}

# run models
epfu.mod <- model_run("epfu") # big brown bat
lano.mod <- model_run("lano") # silver-haried bat
labo.mod <- model_run("labo") # eastern red bat

# Look at estimates
# Overall VIP
epfu.mod$mean.vip
lano.mod$mean.vip
labo.mod$mean.vip

# Median model coefficient and 95% confidence intervals
epfu.mod$o95
lano.mod$o95
labo.mod$o95

# Plot to explore results
m=matrix(1:3,nrow=3,ncol=1)
layout(m)
par(mar=c(5,7,2,2))
par(ps=16)
labels=c("REGIONAL","TREE100","TREE500","TREE1000","OPENVEG100","OPENVEG500","OPENVEG1000","IMP100","IMP500","IMP1000")

bat_plot <- function(species,confint,parpob){
  # big brown bat
  plot(1:10, col="white", xlim=c(-2,5), xlab="", yaxt='n', xaxt='n', ylab="", bty='n', main=species)
  text(-2.5,c(1:10),labels=labels,xpd=TRUE)
  axis(1,at=seq(-2,5,1), labels=FALSE)
  abline(v = 0, lty = 2)
  
  # Points
  for (i in 1:nrow(confint)){
    points(confint[i,2],i,pch=16, col="darkgrey",cex=1.6)
  }
  
  # CI's
  for(i in 1:nrow(confint)){
    segments(confint[i,1],i,confint[i,3], lwd=1.5,col="darkgrey")
    segments(confint[i,1],i,confint[i,3], lwd=3, col="darkgrey")
  }
  
  # Variable Probabilities
  text(4.8,10.75,labels="VIP",font=2,xpd=TRUE)
  for(i in 1:length(parpob)){
    text(4.8,i,labels=parpob[i],xpd=TRUE)
  }
}

# Plot in 3x1 panel
bat_plot("Big brown bat",epfu.mod$o95,epfu.mod$parprob)
bat_plot("Silver-haired bat",lano.mod$o95,lano.mod$parprob)
bat_plot("Eastern red bat",labo.mod$o95,labo.mod$parprob)

# Testing model fit using leave-one-out & McFadden's Psuedo R-square
# First we predict out of sample data point using the final model and a null model
# We then calculate the likelihood that the model can predict that point
# And finally use those likelihoods to calculate McFadden's psuedo R-square values

loo <- function(species, pi) {
  # Choose the detection data for the species to be analyzed (epfu, labo, lano)
  choice <- sp_dat[[species]]$spcount
  
  # Create new data list for this part of the analysis
  # Check where we have a data point to remove
  ch <- which(is.na(choice)==FALSE)
  
  # Create vector to store likelihood of each datapoint
  lik_full <- rep(0, length(ch))
  lik_null <- rep(0, length(ch))
  
  # Create initial values with random starting values for EVERY parameter
  inits_fnc <- function(chain){
    gen_list <- function(chain = chain){
      list( 
        z = z2,
        laplace = runif(ncov, -3, 3),
        int_mu = runif(1, -3,3),
        int_yr_mu = runif(3, -3, 3),
        lp_mu = runif(1, -3, 3),
        lp_mu_yr = runif(3, -3, 3),
        int_tau_session = rgamma(1,1,1),
        int_tau_year = rgamma(1,1,1),
        lp_tau_session = rgamma(1,1,1),
        lp_tau_year = rgamma(1,1,1),
        lp_0 = runif(8, -3, 3),
        beta_0 = runif(8, -3, 3),
        lambda = rgamma(1,1,1),
        pp=rbeta(1,1,1),
        
        .RNG.name = switch(chain,
                           "1" = "base::Wichmann-Hill",
                           "2" = "base::Marsaglia-Multicarry",
                           "3" = "base::Super-Duper",
                           "4" = "base::Mersenne-Twister",
                           "5" = "base::Wichmann-Hill",
                           "6" = "base::Marsaglia-Multicarry",
                           "7" = "base::Super-Duper",
                           "8" = "base::Mersenne-Twister"),
        .RNG.seed = sample(1:1e+06, 1)
      )
    }
    
    return(switch(chain,           
                  "1" = gen_list(chain),
                  "2" = gen_list(chain),
                  "3" = gen_list(chain),
                  "4" = gen_list(chain),
                  "5" = gen_list(chain),
                  "6" = gen_list(chain),
                  "7" = gen_list(chain),
                  "8" = gen_list(chain)
    )
    )
  }
  
  # Loop through each datapoint
  for(i in 1:length(ch)){
    # New dataset to remove one data point
    choice2 <- choice
    
    # Make datapoint NA (it will then be estimated)
    choice2[ch[i]] <- NA
    data_list <- list(y = choice2, x = x,
                      yrs = indx$yrs, sites = indx$sites, 
                      yr_mu_vec =yrmu_vec,
                      ncov = ncov, nyr = 3, nsession = length(yrmu_vec),
                      n = length(choice), jmat =  jmatlong$value,
                      pred_y = choice[ch[i]], loo = ch[i], 
                      pi=pi) # pi must be changed to indicate the final model for each species
    
    # Create starting values for z parameter
    z2 <- choice2
    z2[z2>1] <- 1
    z2[is.na(z2)] <- 0
    
    # Run JAGS model for full model
    m_out <- run.jags(model = "gallo_et_al_2018_base_lasso_loocv_model.R",
                      monitor = c("lik"),
                      inits = inits_fnc,
                      data = data_list,
                      n.chains = 7, adapt = 5000, burnin = 20000, sample = 20000,
                      summarise = FALSE,
                      plots = FALSE,
                      method = "parallel")
    
    # Store likelihood values for each data point
    lik_full[i] <- mean(as.matrix(as.mcmc.list(m_out), chains = TRUE)[,2])
    
    # Run JAGS model for null model
    m_out_null <- run.jags(model = "gallo_et_al_2018_base_lasso_loocv_null_model.R",
                      monitor = c("lik"),
                      inits = inits_fnc,
                      data = data_list,
                      n.chains = 7, adapt = 5000, burnin = 20000, sample = 20000,
                      summarise = FALSE,
                      plots = FALSE,
                      method = "parallel")
    lik_null[i] <- mean(as.matrix(as.mcmc.list(m_out_null), chains = TRUE)[,2])
  }
  
  #Psuedo R-square
  l <- sum(log(lik_full))
  l0 <- sum(log(lik_null))
  r2 <- 1 - (l/l0)
  
  return(list(lik_full=lik_full, mod_full=m_out, lik_null=lik_null, mod_null=m_out_null, r2=r2))
}

# Must enter species and pi is a vectors of 1's and 0's that indicate which parameters to include in final model
# Note order of pi is in the same order as object 'x'. It is not in the same order as previously graphed
# epfu = c(0,0,0,1,0,0,1,0,0,1)
# lano = c(1,0,0,1,0,1,1,1,1,0)
# labo = c(1,0,0,1,0,1,0,1,0,0)

epfu_loo_full <- loo("epfu", c(0,0,0,1,0,0,1,0,0,1))
lano_loo_full <- loo("lano", c(1,0,0,1,0,1,1,1,1,0))
labo_loo_full <- loo("labo", c(1,0,0,1,0,1,0,1,0,0))

# Check R-square value
epfu_loo_full$r2
lano_loo_full$r2
labo_loo_full$r2

# Using our final models for EPFU and LANO we predict across the greater landscape
# Function to extract land cover data at each scale from a grid of points spaced 100m apart across the greater study region

extract_bufferPredict <- function (buff) {
  # load fine scale cover data
  rm <- raster("2010_HRLC_North_Central_Merge.tif")
  
  # load matrix of points
  points <- readOGR(dsn=".",layer="gallo_et_al_2018_grid_100m_reduced")
  points$id <- as.character(points$id)
  
  # setup parallel backend to use many processors
  cl <- makeCluster(7) #set number of cores
  registerDoParallel(cl) # register backend
  
  # loop through each point
  dataLong <- foreach(i=1:length(points),.combine=rbind, .packages=c("raster", "rgdal", "stringr")) %dopar% {
    # extract land cover data for each point, given buffer size
    Landcover <- prop.table(table(extract(rm, points[i,], buffer=buff)))
    if(length(Landcover)==0) {
      Landcover <- NA
      names(Landcover) <- "BLANK"
    }
    # summarize each site's data by proportion of each cover type
    # convert to data frame
    data.frame(id = points[i,]$id,
               cover = names(Landcover),
               percent = as.numeric(Landcover)
    )
  }
  # stop cluster
  stopCluster(cl)
  
  # reshape data
  mydf_reshape <- reshape(dataLong,idvar="id",timevar="cover", direction="wide")
  mydf_reshape$id <- as.character(mydf_reshape$id)
  
  # remove dataLong to save memory
  rm(dataLong)
  
  # NA's to 0 
  mydf_reshape[is.na(mydf_reshape)] <- 0
  mydf_reshape <- mydf_reshape[,-grep("BLANK",colnames(mydf_reshape))]
  
  # create cover name column
  # order by site and order columns 1-7
  df <- mydf_reshape[order((as.numeric(mydf_reshape$id))),order(names(mydf_reshape))]
  
  # setup parallel backend to use many processors
  cl <- makeCluster(7) #set number of cores
  
  # 5,6,& 7 combined are impervious cover
  df$imp <- 0
  df$imp <- parApply(cl,df[,6:8],1,sum)
  
  # stop cluster
  stopCluster(cl)
  
  # remove sites column and proportions of ID #3,6 & 7 because we dont want bare ground, roads and paved areas
  df <- df[,c(1:3,9)]
  colnames(df) <- c("id",paste("tree_", buff, sep=""),paste("openveg_", buff, sep=""), paste("imp_", buff, sep=""))
  
  return(df)
  
}

# Combine objects to make one dataframe (takes several days)
data.comb <- cbind(extract_bufferPredict(100), extract_bufferPredict(500)[,2:4], extract_bufferPredict(1000)[,2:4])

# Inidcate urban vs. not urban points
# Points >50 km from urban center were considered urban
urban.points <- readOGR(dsn=".",layer="gallo_et_al_2018_urban_points")
urban <- as.numeric(as.character(urban.points$id))
data.comb$trt <- 0

# In original grid of points indicate urban sites with a 1
data.comb[which(data.comb$id %in% urban),11] <- 1
# Check that its all set up correctly
table(data.comb$trt)
head(data.comb)

# Scale new covariates by the mean and sd of original covariates
tree.100 <- (data.comb$tree_100 - mean(orig_covs$tree_100))/sd(orig_covs$tree_100)
tree.500 <- (data.comb$tree_500 - mean(orig_covs$tree_500))/sd(orig_covs$tree_500)
tree.1000 <- (data.comb$tree_1000 - mean(orig_covs$tree_1000))/sd(orig_covs$tree_1000)
openveg.100 <- (data.comb$openveg_100 - mean(orig_covs$openveg_100))/sd(orig_covs$openveg_100)
openveg.500 <- (data.comb$openveg_500 - mean(orig_covs$openveg_500))/sd(orig_covs$openveg_500)
openveg.1000 <- (data.comb$openveg_1000 - mean(orig_covs$openveg_1000))/sd(orig_covs$openveg_1000)
imp.100 <- (data.comb$imp_100 - mean(orig_covs$imp_100))/sd(orig_covs$imp_100)
imp.500 <- (data.comb$imp_500 - mean(orig_covs$imp_500))/sd(orig_covs$imp_500)
imp.1000 <- (data.comb$imp_1000 - mean(orig_covs$imp_1000))/sd(orig_covs$imp_1000)
urb_sites <- data.comb$trt

# Function to use later for backtransformation
ex <- function(x) {
  1/(1 + exp(-x))
}

# Function to extract data, predict p, and create raster
get_raster <- function (model_mat) {
  
  md <- model_mat
  
  # calculate variable probabilities
  best_pi <- md[, grep("pi\\[", colnames(md))]
  par_prob <- apply(best_pi, 2, mean)
  names(par_prob) <- colnames(x)
  par_prob <- t(t(par_prob))
  rownames(par_prob) <- labels
  
  # selecting parameters to include in the model
  # relative model selection, if the prob of the param being selected is < the mean prob of all params its centered on zero
  tozero <- which(par_prob<mean(par_prob))
  
  # get 95 % credible intervals
  best <- md[, grep("beta\\[", colnames(md))]
  
  for(i in 1:ncol(best)){
    zero <- which(best[,i]==0)
    if(i %in% tozero){
      best[-zero,i] <- best[-zero,i] - median(best[-zero,i])
    }
  }
  
  # loop to grab values from mcmc chains
  o95 <- matrix(0, ncol = 3, nrow = c(ncol(best)))
  m1 <- rep(0, ncol(best))
  for(i in 1:ncol(best)){
    n0 <- which(best[,i] != 0)
    m1[i] <- length(n0)
    o95[i,] <- quantile(best[n0,i], probs = c(0.025, 0.5, 0.975))
  }
  rownames(o95) <- labels
  
  # get intercept value
  intercept <- md[, grep("lp_mu", colnames(md))]
  b0 <- median(intercept[,1])
  
  # run predictive model
  p.logit <- b0 + o95[1,2]*imp.100 + o95[2,2]*imp.1000 + o95[3,2]*imp.500 + o95[4,2]*openveg.100 + o95[5,2]*openveg.1000 + o95[6,2]*openveg.500 +
    + o95[7,2]*tree.100 + o95[8,2]*tree.1000 + o95[9,2]*tree.500 + o95[10,2]*urb_sites
  
  # transform to probability scale
  p <- ex(p.logit)
  
  # matrix dimensions
  n.col <- 852 #number of points in first row
  a <- length(urb_sites) # total points (area)
  n.row <- a/n.col
  
  # fill matrix
  mat <- matrix(p,nrow=n.row,ncol=n.col, byrow=TRUE)
  
  # create raster
  crs <- CRS("+init=epsg:26916")
  p.raster <- raster(mat, crs=crs)
  extent(p.raster) <- c(372021.50,457121.50,4619577.23,4679077.23)
  
  return (list(b0=b0,o95=o95,raster=p.raster,mat=mat))
  
}

# extract info
epfu.info <- get_raster(epfu.mod$model_mat)
labo.info <- get_raster(labo.mod$model_mat)

# Plot and map predictive rasters

# Features to add to maps
chicago.river <- readOGR(dsn=".",layer="ChicagoRiver_clip2")
fox.river <- readOGR(dsn=".",layer="FoxRiver_Clip")
lake <- readOGR(dsn=".",layer="lake_michigan_shoreline")

# Plot rasters
m <- matrix(1:3,nrow=1,ncol=3,byrow=TRUE)
bw.colramp <- c('#000000','#1c1c1c','#343434','#4b4b4b','#666666','#808080','#9c9c9c','#b8b8b8','#d6d6d6','#f5f5f5')
my.colramp <- c('#0000ff','#4565ff','#6297fe','#77b9fc','#85ccfa','#f6d57c','#feb15c','#ff8c3e','#ff611f','#ff0000')
tiff("prediction_plot.tiff", height = 6, width = 28, units = "in", res = 600, compression = "lzw")
layout(m)
par(mar=c(5,4,4,12) + 0.1)
par(ps=22)

# EPFU
brks=seq(0,1,0.1)
plot(epfu.info$raster, breaks=brks,col=my.colramp,box=FALSE,axes=FALSE, legend=FALSE, yaxt='n',xaxt='n', xlim=c(380000,460000),ylim=c(4620000,4680000), main="Big brown bat")
plot(lake, col="white",border=FALSE, add=TRUE)
plot(chicago.river,lwd=2,add=TRUE)
plot(fox.river,lwd=2,add=TRUE)
rect(xleft=360000,ybottom=4620000,xright=372000,ytop=4680000, col="white",border=NA)
axis(2,at=c(4620000,4640000,4660000,4680000),line=-2,xpd=TRUE)
axis(1,at=c(380000,400000,420000,440000,460000),line=1)

# LABO
plot(labo.info$raster, breaks=brks,col=my.colramp,box=FALSE,axes=FALSE,yaxt='n',xaxt='n', legend=FALSE,xlim=c(380000,460000),ylim=c(4620000,4680000), main="Eastern red bat")
plot(lake, col="white",border=FALSE, add=TRUE)
plot(chicago.river,lwd=2,add=TRUE)
plot(fox.river,lwd=2,add=TRUE)
rect(xleft=360000,ybottom=4620000,xright=372000,ytop=4680000, col="white",border=NA)
axis(1,at=c(380000,400000,420000,440000,460000),line=1)

# Combine the two species into a prob raster
sp.mat <- epfu.info$mat*labo.info$mat
crs <- CRS("+init=epsg:26916")
com.raster <- raster(sp.mat, crs=crs)
extent(com.raster) <- c(372021.50,457121.50,4619577.23,4679077.23)
# Plot
plot(com.raster, breaks=brks,col=my.colramp,box=FALSE,axes=FALSE,yaxt='n',xaxt='n', xlim=c(380000,460000),ylim=c(4620000,4680000), legend.args=list(text='Pr(activity)', side=2, font=2, line=2.5, cex=0.8),main="Species combined")
plot(lake, col="white",border=FALSE, add=TRUE)
plot(chicago.river,lwd=2,add=TRUE)
plot(fox.river,lwd=2,add=TRUE)
rect(xleft=360000,ybottom=4620000,xright=372000,ytop=4680000, col="white",border=NA)
axis(1,at=c(380000,400000,420000,440000,460000),line=1)

dev.off()

# Caculate summary stats
# Proportion of study area with high and low activity
prop.table(table(epfu.info$mat>0.69))
prop.table(table(epfu.info$mat<0.31))

prop.table(table(labo.info$mat>0.69))
prop.table(table(labo.info$mat<0.31))

# Mean activity across study area
mean(epfu.info$mat)
mean(labo.info$mat)
mean(sp.mat)

# Calculate summary statistics for data collection
bats <- readRDS("gallo_et_al_2018_bat_detection_FullData.RDS")
colnames(bats) <- tolower(colnames(bats))

# Total number of species
unique(bats$decision) 

# Number of species in urban and rural sites
urban.sites <- covs[which(covs$trt==1),1] # urban
rural.sites <- covs[which(covs$trt==0),1] # rural
unique(bats$decision[which(bats$site %in% urban.sites & bats$spcount!=0)]) # urban
unique(bats$decision[which(bats$site %in% rural.sites& bats$spcount!=0)]) # rural

# Most common species
summary.table <- cbind(rbind(sum(bats$spcount[which(bats$decision=="Epfu" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Labo" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Laci" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Lano" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Mylu" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Myse" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Nyhu" & bats$site %in% urban.sites)]),
                            sum(bats$spcount[which(bats$decision=="Pesu" & bats$site %in% urban.sites)])),
                      rbind(sum(bats$spcount[which(bats$decision=="Epfu" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Labo" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Laci" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Lano" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Mylu" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Myse" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Nyhu" & bats$site %in% rural.sites)]),
                            sum(bats$spcount[which(bats$decision=="Pesu" & bats$site %in% rural.sites)]))
)

colnames(summary.table) <- c("Urban","Rural")
rownames(summary.table) <- unique(bats$decision)

# Combine myotis species
rownames(summary.table)[5] <- "Myotis"
summary.table[5,1] <- summary.table[5,1] + 1
summary.table <- summary.table[-6,]

