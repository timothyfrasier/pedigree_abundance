##################################
#        abundanceSIM            #
#                                #
# This function takes input      #
# parameters and simulates a     #
# population with a defined      #
# number of mothers, fathers,    #
# and offspring. It then samples #
# those and returns data on the  #
# characteristics of those       #
# sampled individuals, as needed #
# for abundance inference.       #
#--------------------------------#
# REQUIRES:                      #
#  N = population size           #
#                                #
#  propReproductive = proportion #
#    of population that is       #
#    reproductive.               #
#                                #
#  propSampled = proportion of   #
#    population that is sampled. #
#--------------------------------#
# RETURNS:                       #
#   Ns = Number of sampled       #
#        individuals.            #
#   Nin = Number of inferred     #
#         individuals.           #
#   Bs = Number of sampled       #
#        breeders.               #
##################################

abundanceSim = function(N, propReproductive, propSampled) {

  #--- Create Reference Population ---#
  nOffspring = round((N * propReproductive) / 2)
  offspring = 1:nOffspring
  moms = (nOffspring + 1):(nOffspring * 2)
  dads = ((nOffspring * 2) + 1):(nOffspring * 3)
  others = ((nOffspring * 3) + 1):N

  #--- Create and Sample Population ---#
  population = 1:N
  population = sample(population, size = length(population), replace = FALSE)

  sampSize = round(propSampled * N)
  sampled = population[1:sampSize]

  #--- Collect Data ---#
  Ns = sampSize
  momSampled = moms %in% sampled
  dadSampled = dads %in% sampled
  offSampled = offspring %in% sampled

  triads = 0
  mDyads = 0
  dDyads = 0
  justOff = 0

  for (i in 1:nOffspring) {
    if (offSampled[i] == TRUE & momSampled[i] == TRUE & dadSampled[i] == TRUE) {
      triads = triads + 1
   } else if (offSampled[i] == TRUE & momSampled[i] == TRUE & dadSampled[i] == FALSE) {
      mDyads = mDyads + 1
    } else if (offSampled[i] == TRUE & momSampled[i] == FALSE & dadSampled[i] == TRUE) {
      dDyads = dDyads + 1
    } else if (offSampled[i] == TRUE & momSampled[i] == FALSE & dadSampled[i] == FALSE) {
      justOff = justOff + 1
    }
  }

  Bs = (triads * 2) + mDyads + dDyads
  Nin = mDyads + dDyads + (2 * justOff)
  
  results = c(Ns, Nin, Bs)

  return(results)
}  
  