# Read the command line arguments
args <- commandArgs(TRUE)
Env  <- args[1]
sim  <- args[2]
rep  <- args[3]

# load the user's environment
load(Env)

# Figure out the appropriate output file location
Output.File <- paste(..Sim.Directory, '/OutputFiles/sim',sim,'rep',rep,'.RData', sep='')

if( !file.exists(Output.File) ){
  # if there were any local libraries the user had, include them
  # in the library search path.
  .libPaths( ..LibPaths )

  # load the packages that were loaded when the user called us
  lapply(..LoadedPackages, require, character.only = TRUE)

  sim <- do.call(what = ..Sim.Function, 
                 args = lapply(..Params[sim,], identity))
  save('sim', file=Output.File)
}
                            