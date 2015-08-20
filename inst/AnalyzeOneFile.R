# 1) Loads the R envirnment in the input file 
# 2) calls the user function (in the loaded envirnment)
# 3) saves the output as "sim"
args <- commandArgs(TRUE)
Input.File  <- args[1]
load(Input.File)

# if there were any local libraries the user had, include them
# in the library search path.
.libPaths( .LibPaths )

# load the packages that were loaded when the user called us
lapply(..LoadedPackages, require, character.only = TRUE)

if( !file.exists(Output.File) ){
   sim <- do.call(what = ..Sim.Function, 
                  args = lapply(..Params, identity))
   save('sim', file=Output.File)
}
                            