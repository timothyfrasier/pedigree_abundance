####################################
#     pedigreeAbundance.R          #
#                                  #
# This code estimates abundance    #
# via pedigree reconstruction      #
# using a method modified from     #
# Creel & Rosenblatt (2013)        #
# Ecology and Evolution 3(5):      #
# 1294-1304.                       #
#                                  #
# Requires:                        #
#    Ns = # sampled                #
#    Nin = # inferred via pedigree #
#    Bs = # of breeders sampled    #
####################################


pedigreeAbundance = function(Ns, Nin, Bs, mode = "singular") {
  
  #------------------------------#
  # Load Approapriate libraries  #
  # and things.                  #
  #------------------------------#
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  #----------------------------------------#
  # Probability of Being Sampled.          #
  #----------------------------------------#
  Psampled = rep(0, times = (Bs + Nin))
  
  for (i in 1:Bs) {
    Psampled[i] = 1
  }
  Psampled = sample(Psampled, size = length(Psampled), replace = FALSE)
  N = length(Psampled)
  
  #-------------------------------#
  # Write data as a list for STAN #
  #-------------------------------#
  dataList = list (
    y = Psampled,
    N = N
  )
  
  #------------------------------#
  #       Define the Model       #
  #------------------------------#
  modelString = "
    data {
      int<lower=0> N;
      int y[N];
    }

    parameters {
      real theta;
    }

    model {
      // Likelihood
      y ~ bernoulli(theta);

      // Priors
      theta ~ uniform(0, 1);
    }
    "
  writeLines(modelString, con="model.stan")
  
  #--------------------------#
  #     run STAN             #
  #--------------------------#
  stanFit = stan(file = "model.stan", 
                 data = dataList, 
                 pars = "theta",
                 warmup = 2000,
                 iter = 10000, 
                 chains = 3)
  
  #-------------------------------#
  # Extract & Summarize the Data  #
  #-------------------------------#
  mcmcChains = as.data.frame(stanFit)
  Ps = mcmcChains[, "theta"]
  Pns = 1 - Ps
  modePns = getmode(Pns)
  lowPns = as.numeric(quantile(Pns, probs = 0.025))
  highPns = as.numeric(quantile(Pns, probs = 0.975))
  
  #----------------------------#
  #   Plot Data If Desired     #
  #----------------------------#
  if (mode == "singular") {
    histInfo = plotPost(Ps, xlab = expression(paste(P[sampled])), showMode = TRUE)
    histInfo = plotPost(Pns, xlab = expression(paste(P[not_sampled])), showMode = TRUE)
  }
  
  #----------------------------------------#
  # Make prob. of not breeding (Pnb) a     # 
  # distribution and estimate with model.  #
  # We actually estimate P(breeding) with  #
  # the model and then subtract that from  #
  # 1 to get P(not breeding).              #
  #----------------------------------------#
  Pbreeding = rep(0, times = Ns)
  
  for (i in 1:(Bs + Nin)) {
    Pbreeding[i] = 1
  }
  Pbreeding = sample(Pbreeding, size = length(Pbreeding), replace = FALSE)
  N = length(Pbreeding)
  
  #-------------------------------#
  # Write data as a list for STAN #
  #-------------------------------#
  dataList = list (
    y = Pbreeding,
    N = N
  )
  
  #--------------------------#
  #     run STAN             #
  #--------------------------#
  stanFit = stan(file = "model.stan", 
                 data = dataList, 
                 pars = "theta",
                 warmup = 2000,
                 iter = 10000, 
                 chains = 3)
  
  #-------------------------------#
  # Extract & Summarize the Data  #
  #-------------------------------#
  mcmcChains = as.data.frame(stanFit)
  Pbd = mcmcChains[, "theta"]
  Pnb = 1 - Pbd
  modePbd = getmode(Pbd)
  lowPbd = as.numeric(quantile(Pbd, probs = 0.025))
  highPbd = as.numeric(quantile(Pbd, probs = 0.975))
  
  #----------------------------#
  #   Plot Data If Desired     #
  #----------------------------#
  if (mode == "singular") {
    histInfo = plotPost(Pbd, xlab = expression(paste(P[bd])), showMode = TRUE)
    histInfo = plotPost(Pnb, xlab = expression(paste(P[nb])), showMode = TRUE)
  }
  
  
  #-------------------------------#
  # Non-sampled non-breeders      #
  #-------------------------------#
  Nnot_sampled = round((Ns / Ps) - Ns)
  modeNnot_sampled = getmode(Nnot_sampled)
  lowNnot_sampled = as.numeric(quantile(Nnot_sampled, probs = 0.025))
  highNnot_sampled = as.numeric(quantile(Nnot_sampled, probs = 0.975))
  
  Nnsnb = round(Nnot_sampled * Pnb)
  Nnsnb = ifelse(Nnsnb < 0, 0, Nnsnb)
  modeNnsnb = getmode(Nnsnb)
  lowNnsnb = as.numeric(quantile(Nnsnb, probs = 0.025))
  highNnsnb = as.numeric(quantile(Nnsnb, probs = 0.975))
  
  #----------------------------#
  #   Plot Data If Desired     #
  #----------------------------#
  if (mode == "singular") {
    histInfo = plotPost(Nnot_sampled, xlab = expression(paste(N[not_sampled])), showMode = TRUE)
    histInfo = plotPost(Nnsnb, xlab = expression(paste(N[nsnb])), showMode = TRUE)
  }
  
  
  #--------------------------------------#
  # Non-sampled, non-inferred, breeders  #
  #--------------------------------------#
  Nbns = round(Nnot_sampled * Pbd)
  Nbns = ifelse(Nbns < 0, 0, Nbns)
  modeNbns = getmode(Nbns)
  lowNbns = as.numeric(quantile(Nbns, probs = 0.025))
  highNbns = as.numeric(quantile(Nbns, probs = 0.975))
  
  Nbnsni = round(Nbns - Nin)
  Nbnsni = ifelse(Nbnsni < 0, 0, Nbnsni)
  modeNbnsni = getmode(Nbnsni)
  lowNbnsni = as.numeric(quantile(Nbnsni, probs = 0.025))
  highNbnsni = as.numeric(quantile(Nbnsni, probs = 0.975))
  
  #----------------------------#
  #   Plot Data If Desired     #
  #----------------------------#
  if (mode == "singular") {
    histInfo = plotPost(Nbns, xlab = expression(paste(N[bns])), showMode = TRUE)
    histInfo = plotPost(Nbnsni, xlab = expression(paste(N[bnsni])), showMode = TRUE)
  }
  
  
  #------------------#
  #   Estimate N     #
  #------------------#
  N = Ns + Nin + Nnsnb + Nbnsni
  modeN = getmode(N)
  lowN = as.numeric(quantile(N, probs = 0.025))
  highN = as.numeric(quantile(N, probs = 0.975))
  
  #----------------------------#
  #   Plot Data If Desired     #
  #----------------------------#
  if (mode == "singular") {
    histInfo = plotPost(N, xlab = "N", showMode = TRUE)
  }
  
  if (mode == "singular") {
    results = list(Pns, Pbd, Nnot_sampled, Nnsnb, Nbns, Nbnsni, N)
  } else {
    results = c(modePns, lowPns, highPns, modePbd, lowPbd, highPbd, modeNnot_sampled, lowNnot_sampled, highNnot_sampled, modeNnsnb, lowNnsnb, highNnsnb, modeNbns, lowNbns, highNbns, modeNbnsni, lowNbnsni, highNbnsni, modeN, lowN, highN)
  }  
  
  return(results)
}


#####################################
#    pedigreeAbundance_h            #
#-----------------------------------#
# This function estimates abundance #
# for MULTIPLE years/time periods/  #
# locations simultaneously using a  #
# hierarchical Bayesian model via   #
# pedigree reconstruction using a   #
# method modified from Creel &      #
# Rosenblatt (2013) Ecology and     #
# Evolution 3(5): 1294-1304.        #
#                                   #
# Requires:                         #
#    Ns = a vector of the # sampled #
#         individuals for each      #
#         "event" (time/location)   #
#    Nin = a vector of the # of     # 
#          inferred individuals for #
#          each "event" via pedigree#
#          reconstruction.          #
#    Bs = a vector of the # of      #
#         breeders sampled for each #
#         "event"                   #
#####################################
pedigreeAbundance_h = function (Ns, Nin, Bs) {

  #----------------------------------#
  # Load libraries and organize data #
  #----------------------------------#
  library(rstan)
  options(mc.cores = parallel::detectCores())
  library(dplyr)
  
  
  #----------------------------------#
  # Check the data for same length   #
  #----------------------------------#
  if (length(Ns) != length(Nin) | length(Ns) != length(Bs) | length(Nin) != length(Bs)) {
    stop("ERROR: Must have same number of data points for Ns, Nin, and Bs")
  }
  

  #--------------#
  #   Psampled   #
  #--------------#
  total = Bs + Nin
  success = Bs
  events = length(total)
  
  percentage = success / total
  PsMean = mean(percentage)
  PsSD = sd(percentage)


  #-------------------------------#
  # Write data as a list for STAN #
  #-------------------------------#
  dataList = list (
    total = total,
    success = success,
    events = events,
    PsMean = PsMean,
    PsSD = (PsSD * 3)
  )


  #------------------------------#
  #       Define the Model       #
  #------------------------------#
  modelString = "
    data {
      int<lower=0> events;
      int total[events];
      int success[events];
      real PsMean;
      real PsSD;
    }

    parameters {
      real<lower=0, upper=1> theta[events];
      real<lower=0, upper=1> eventMean;
    }

    model {
      // Likelihood
      for (i in 1:events) {
        success[i] ~ binomial(total[i], theta[i]);
      }  

      // Priors
      for (i in 1:events) {
        theta ~ cauchy(eventMean, PsSD); 
      }
    
      // Hyperpriors
      eventMean ~ normal(PsMean, PsSD);
    }
  "
  writeLines(modelString, con="model.stan")

  #--------------------------#
  #     run STAN             #
  #--------------------------#
  stanFit = stan(file = "model.stan", 
                 data = dataList, 
                 pars = c("theta", "eventMean"),
                 warmup = 2000,
                 iter = 10000, 
                 chains = 3)

  print(stanFit)
  mcmcChains = as.data.frame(stanFit)

  Psampled = select(mcmcChains, -(eventMean))
  Ps = select(Psampled, -(lp__))
  Pns = 1 - Ps
  meanPs = as.numeric(apply(Ps, 2, mean))
  lowPs = as.numeric(apply(Ps, 2, quantile, probs = 0.025))
  highPs = as.numeric(apply(Ps, 2, quantile, probs = 0.975))


  #----------------------------------------#
  # Make prob. of not breeding (Pnb) a     # 
  # distribution and estimate with model.  #
  # We actually estimate P(breeding) with  #
  # the model and then subtract that from  #
  # 1 to get P(not breeding).              #
  #----------------------------------------#
  total = Ns
  success = Bs + Nin
  events = length(total)

  percentage = success / total
  PbMean = mean(percentage)
  PbSD = sd(percentage)

  
  #-------------------------------#
  # Write data as a list for STAN #
  #-------------------------------#
  dataList = list (
    total = total,
    success = success,
    events = events,
    PbMean = PbMean,
    PbSD = (PbSD * 3)
  )

  #------------------------------#
  #       Define the Model       #
  #------------------------------#
  modelString = "
    data {
      int<lower=0> events;
      int total[events];
      int success[events];
      real PbMean;
      real PbSD;
    }

    parameters {
      real<lower=0, upper=1> theta[events];
      real<lower=0, upper=1> eventMean;
    }

    model {
      // Likelihood
      for (i in 1:events) {
        success[i] ~ binomial(total[i], theta[i]);
      }  

      // Priors
      for (i in 1:events) {
        theta ~ cauchy(eventMean, PbSD); 
      }
    
      // Hyperpriors
      eventMean ~ normal(PbMean, PbSD);
    }
  "
  writeLines(modelString, con="model.stan")

  #--------------------------#
  #     run STAN             #
  #--------------------------#
  stanFit = stan(file = "model.stan", 
                 data = dataList, 
                 pars = c("theta", "eventMean"),
                 warmup = 2000,
                 iter = 10000, 
                 chains = 3)

  print(stanFit)
  mcmcChains = as.data.frame(stanFit)

  Pbreeding = select(mcmcChains, -(eventMean))
  Pbd = select(Pbreeding, -(lp__))
  Pnb = 1 - Pbd
  meanPbd = as.numeric(apply(Pbd, 2, mean))
  lowPbd = as.numeric(apply(Pbd, 2, quantile, probs = 0.025))
  highPbd = as.numeric(apply(Pbd, 2, quantile, probs = 0.975))


  #-------------------------------#
  # Non-sampled non-breeders      #
  #-------------------------------#
  Nnot_sampled = data.frame(matrix(NA, nrow = length(Ps[, 1]), ncol = events))
  for (i in 1:length(Ns)) {
    for (j in 1:length(Ps[, 1])) {
      Nnot_sampled[j, i] = round((Ns[i] / Ps[j, i]) - Ns[i])
    }
  }

  meanNnot_sampled = as.numeric(apply(Nnot_sampled, 2, mean))
  lowNnot_sampled = as.numeric(apply(Nnot_sampled, 2, quantile, probs = 0.025))
  highNnot_sampled = as.numeric(apply(Nnot_sampled, 2, quantile, probs = 0.975))

  Nnsnb = round(Nnot_sampled * Pnb)
  for (i in 1:length(Nnsnb[, 1])) {
    for (j in 1:length(Nnsnb[1, ])) {
      if (Nnsnb[i, j] < 0) {
        Nnsnb[i, j] = 0
      }
    }
  }
  
  meanNnsnb = as.numeric(apply(Nnsnb, 2, mean))
  lowNnsnb = as.numeric(apply(Nnsnb, 2, quantile, probs = 0.025))
  highNnsnb = as.numeric(apply(Nnsnb, 2, quantile, probs = 0.975))


  #--------------------------------------#
  # Non-sampled, non-inferred, breeders  #
  #--------------------------------------#
  Nbns = round(Nnot_sampled * Pbd)
  for (i in 1:length(Nbns[, 1])) {
    for (j in 1:length(Nbns[1, ])) {
      if (Nbns[i, j] < 0) {
        Nbns[i, j] = 0
      }
    }
  }

  meanNbns = as.numeric(apply(Nbns, 2, mean))
  lowNbns = as.numeric(apply(Nbns, 2, quantile, probs = 0.025))
  highNbns = as.numeric(apply(Nbns, 2, quantile, probs = 0.975))
    
  Nbnsni = data.frame(matrix(NA, nrow = length(Ps[, 1]), ncol = events))
  for (i in 1:length(Nin)) {
    for (j in 1:length(Nbnsni[, 1])) {
      Nbnsni[j, i] = Nbns[j, i] - Nin[i]
    }
  }

  for (i in 1:length(Nbnsni[, 1])) {
    for (j in 1:length(Nbnsni[1, ])) {
      if (Nbnsni[i, j] < 0) {
        Nbnsni[i, j] = 0
      }
    }
  }
  
  meanNbnsni = as.numeric(apply(Nbnsni, 2, mean))
  lowNbnsni = as.numeric(apply(Nbnsni, 2, quantile, probs = 0.025))
  highNbnsni = as.numeric(apply(Nbnsni, 2, quantile, probs = 0.975))
    
  #------------------#
  #   Estimate N     #
  #------------------#
  N = data.frame(matrix(NA, nrow = length(Ps[, 1]), ncol = events))
  for (i in 1:length(Nin)) {
    for (j in 1:length(N[, 1])) {
      N[j, i] = Ns[i] + Nin[i] + Nnsnb[j, i] + Nbnsni[j, i]
    }
  }

  meanN = as.numeric(apply(N, 2, mean))
  lowN = as.numeric(apply(N, 2, quantile, probs = 0.025))
  highN = as.numeric(apply(N, 2, quantile, probs = 0.975))

  #------------------------------#
  #   Combine the Data           #
  #------------------------------#
  results = cbind(meanN, lowN, highN, meanNbnsni, lowNbnsni, highNbnsni, meanNbns, lowNbns, highNbns, meanNnsnb, lowNnsnb, highNnsnb, meanPbd, lowPbd, highPbd, meanPs, lowPs, highPs)
  return(results)
}  


############################
#       getmode            #
#--------------------------#
# FUNCTION TO CALCULATE    #
# MODE IN R.               #
############################
getmode = function(x) {
  uniq = unique(x)
  uniq[which.max(tabulate(match(x, uniq)))]
}


#######################################
#             HDIofMCMC               #
#_____________________________________#
# Function for calculating highest    #
# Density intervals (HDIs) from MCMC  #
# data. Required for plotPost.R       #
# function. From Kruschke (2015)      #
# Doing Bayesian Analysis: A Tutorial #
# With R, JAGS, and Stan.             #
#######################################
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}


#######################################
#             plotPost                #
#_____________________________________#
# Function for calculating highest    #
# Density intervals (HDIs) from MCMC  #
# data. Required for plotPost.R       #
# function. From Kruschke (2015)      #
# Doing Bayesian Analysis: A Tutorial #
# With R, JAGS, and Stan.             #
#######################################
plotPost = function( paramSampleVec , credMass=0.95 , compVal=NULL ,
                     HDItextPlace=0.7 , ROPE=NULL , yaxt=NULL , ylab=NULL ,
                     xlab=NULL , cex.lab=NULL , cex=NULL , xlim=NULL , main=NULL ,
                     col=NULL , border=NULL , showMode=F , showCurve=F , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Parameter"
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"
  
  postSummary = matrix( NA , nrow=1 , ncol=11 , 
                        dimnames=list( c( xlab ) , 
                                       c("mean","median","mode",
                                         "hdiMass","hdiLow","hdiHigh",
                                         "compVal","pcGTcompVal",
                                         "ROPElow","ROPEhigh","pcInROPE")))              
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  if ( is.null(breaks) ) {
    breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                     by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display mean or mode:
  if ( showMode==F ) {
    meanParam = mean( paramSampleVec )
    text( meanParam , cenTendHt ,
          bquote(mean==.(signif(meanParam,5))) , adj=c(.5,0) , cex=cex )
  } else {
    dres = density( paramSampleVec )
    modeParam = dres$x[which.max(dres$y)]
    text( modeParam , cenTendHt ,
          bquote(mode==.(signif(modeParam,5))) , adj=c(.5,0) , cex=cex )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    cvCol = "darkgreen"
    pcgtCompVal = round( 100 * sum( paramSampleVec > compVal )
                         / length( paramSampleVec )  , 1 )
    pcltCompVal = 100 - pcgtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) ,
           lty="dotted" , lwd=1 , col=cvCol )
    text( compVal , cvHt ,
          bquote( .(pcltCompVal)*"% < " *
                    .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ) ,
          adj=c(pcltCompVal/100,0) , cex=0.8*cex , col=cvCol )
    postSummary[,"compVal"] = compVal
    postSummary[,"pcGTcompVal"] = ( sum( paramSampleVec > compVal ) 
                                    / length( paramSampleVec ) )
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    ropeCol = "darkred"
    pcInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                 / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pcInROPE))*"% in ROPE" ) ,
          adj=c(.5,0) , cex=1 , col=ropeCol )
    
    postSummary[,"ROPElow"]=ROPE[1] 
    postSummary[,"ROPEhigh"]=ROPE[2] 
    postSummary[,"pcInROPE"]=pcInROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 )
  text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],5))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],5))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}

