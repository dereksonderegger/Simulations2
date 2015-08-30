# Read the command line arguments
args <- commandArgs(TRUE)
Env  <- args[1]
sim  <- as.integer(args[2])
rep  <- as.integer(args[3])

# load the user's environment
load(Env)

# Figure out the appropriate output file location
Output.File <- paste(..Sim.Directory, '/OutputFiles/sim',sim,'rep',rep,'.RData', sep='')

if( !file.exists(Output.File) ){
  # if there were any local libraries the user had, include them
  # in the library search path.
  .libPaths( ..LibPaths )

  # load the packages that were loaded when the user called us
  suppressPackageStartupMessages(
    lapply(..LoadedPackages, require, character.only = TRUE, quietly=TRUE, warn.conflicts=FALSE)
  )
  
  # set up the random number generator
  #.lec.CurrentStream(paste('sim',sim,'rep',rep, sep=''))
  set.seed(..RNG.seeds[[paste('sim',sim,'rep',rep)]])
  # Below is a truly nasty hack...
  # set.seed(.Random.seed + 10000*sim + rep)
  
  # call the user's function
  sim <- do.call(what = ..Sim.Function, 
                 args = lapply(..Params[sim,], identity))
  save('sim', file=Output.File)
}
                            