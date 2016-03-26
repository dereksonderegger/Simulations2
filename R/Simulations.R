#' @export
prep.sims <- function( sim.function, param.matrix, 
                       num.reps=10, sim.directory='sim_directory',
                       Est.Time='20:00'){

	num.sims <- dim(param.matrix)[1]

	# create the sim_directory if necessary
	if( !file.exists(sim.directory) ){ 
 		dir.create(sim.directory)
	}  
	if( !file.exists(paste(sim.directory,'/OutputFiles',sep='')) ){
	  dir.create(paste(sim.directory,'/OutputFiles',sep=''))
	}
	if( !file.exists(paste(sim.directory,'/ConsoleFiles',sep='')) ){
	  dir.create(paste(sim.directory,'/ConsoleFiles',sep=''))
	}


  # record what is in our environment, packages, and local libraries
	Sim.Env <- globalenv()
	Sim.Env$..LoadedPackages <- (.packages())
	Sim.Env$..LibPaths       <- .libPaths()
	Sim.Env$..Sim.Function   <- sim.function
	Sim.Env$..Sim.Directory  <- sim.directory
	Sim.Env$..Params         <- param.matrix
	Sim.Env$..RNG.seeds      <- make.RNG.seeds(.Random.seed, num.sims, num.reps)
	Sim.Env$..RNG.seed.names <- make.RNG.seed.names(num.sims, num.reps)
	
	save(list  = ls(all.names=TRUE, envir=Sim.Env),
	     file  = paste(sim.directory,'/Env.RData',sep=''),
	     envir = Sim.Env )


	# Generate files that contains one command for rep
	script <- system.file('AnalyzeOneFile.R', package='Simulations2')
	for(i in 1:num.sims){
	  file.create(paste(sim.directory,'/ConsoleFiles/CommandFile_',i,'.txt', sep=''))
	  CommandFile <- file(paste(sim.directory,'/ConsoleFiles/CommandFile_',i,'.txt',sep=''), open='a' )
	  for(j in 1:num.reps){
	    writeLines(str_c(R.home(),'/Rscript ', script, ' Env.RData ', i,' ',j), CommandFile)
	  }
	  close(CommandFile)
	}
	
	# Merge all the sims together into one command file
	#   - If the job is huge, SLURM makes me break it into 1000 job chunks
	#   - If it is small, I want to just have one command file
# 	old.wd <- getwd()
# 	setwd(sim.directory)
# 	call <- paste('cat', paste(paste('CommandFile_', 1:num.sims,'.txt', sep=''), collapse=' '), '> CommandFile.txt')
# 	try(system(call))
#   setwd(old.wd)
  
	# Generate a bash shell file that controls each Sim array
	for( i in 1:num.sims){
  	file.create(paste(sim.directory,'/ConsoleFiles/Driver_',i,'.sh',sep=''))
  	Driver <- file(paste(sim.directory,'/ConsoleFiles/Driver_',i,'.sh',sep=''), open='a' )
  	writeLines('#!/bin/sh', Driver)
  	writeLines(paste('#SBATCH --job-name=Simulation'), Driver)
  	writeLines(paste('#SBATCH --output="',sim.directory,'/ConsoleFiles/Output_%A_%a.txt"',sep=''), Driver)
  	writeLines(paste('#SBATCH --workdir="',sim.directory,'"', sep=''), Driver)
  	writeLines(paste('#SBATCH --array=1-', num.reps, sep=''), Driver)
  	writeLines(paste('#SBATCH --time=', Est.Time, sep=''), Driver)
  	writeLines(paste('module load R', sep=''), Driver)
  	writeLines(paste('command=$(awk "NR==$SLURM_ARRAY_TASK_ID" CommandFile_',i,'.txt)', sep=''), Driver)
  	writeLines(paste('srun $command'), Driver)
	  close(Driver)
	}
# 	file.create(paste(sim.directory,'/Driver.sh',sep=''))
# 	Driver <- file(paste(sim.directory,'/Driver.sh',sep=''), open='a' )
# 	writeLines('#!/bin/sh', Driver)
# 	writeLines(paste('#SBATCH --job-name=Simulation'), Driver)
# 	writeLines(paste('#SBATCH --output="',sim.directory,'/ConsoleFiles/Output_%A_%a.txt"',sep=''), Driver)
# 	writeLines(paste('#SBATCH --workdir="',sim.directory,'"', sep=''), Driver)
# 	writeLines(paste('#SBATCH --array=1-', num.sims*num.reps, sep=''), Driver)
# 	writeLines(paste('#SBATCH --time=', Est.Time, sep=''), Driver)
# 	writeLines(paste('module load R', sep=''), Driver)
# 	writeLines(paste('command=$(awk "NR==$SLURM_ARRAY_TASK_ID" CommandFile.txt)', sep=''), Driver)
# 	writeLines(paste('srun $command'), Driver)
# 	close(Driver)
} 

#' @export
run.sims <- function(sim.directory, SLURM=FALSE){
  old.dir <- getwd()
  setwd(sim.directory)   

  all.Console.Files <- list.files(str_c(sim.directory,'/ConsoleFiles')) 
  all.Command.Files <- all.Console.Files[ str_detect(all.Console.Files, fixed('CommandFile'))]
  all.Driver.Files  <- all.Console.Files[ str_detect(all.Console.Files, fixed('Driver'))]
  
  if(!SLURM){
    all.lines <- NULL
    for(commandfile in all.Command.Files){
      con=file(str_c(sim.directory,'/ConsoleFiles/',commandfile)) 
      lines=readLines(con) 
      close(con)
      all.lines <- c(all.lines, lines)
    }
#    script <- system.file('AnalyzeOneFile.R', package='Simulations2')
    foreach(line = all.lines) %dopar% {
      foo <- system( line ) 
    }
    return(invisible(TRUE))
  }
  else{
    setwd(sim.directory)   
    for( driver in all.Driver.Files){
      system(str_c('sbatch ', driver))
    }
    return(invisible(TRUE))
  }
}

  

#' Summarize a simulation study 
#' 
#' Summarize a simulation study where we have a matrix
#' of parameter values that we are interested. 
#' 
#' @param summary.function A function to summarize a single
#' simulation into a vector of output statistics.
#' @param params A data frame of parameter values. This should
#' be the same data frame as was used to create the simulation.
#' @param sim.directory The local directory to store the simulation
#' results.
#' @examples 
#' # Simulation is creating data and fitting 
#' # a regression model to estimate a slope.
#' Sim.Function <- function(N, alpha, beta, sigma){
#'   x <- runif(N, 0, 10)
#'   y <- alpha + beta*x + rnorm(N, 0, sigma)
#'   model <- lm(y~x)
#'   return(model)
#' }
#' 
#' # Create a matrix where each row represents a set of parameters
#' Params <- expand.grid(N=20, alpha=0, beta=c(0,1), sigma=c(.2, 1))
#' 
#' # Run the simulations.  
#' run.sims(Sim.Function, Params, num.reps=10)
#'
#' # Now that the simulations are run, we want to analyze
#' # them.  Create a function to be applied to each simulation
#' # that results in a single row.  Those rows will be concatenated
#' # together into a large data frame that summarizes the simulation.
#' 
#' # Because the output of the Sim.Function was a linear model output,
#' # that is the input for this function
#' Summary.Function <- function(model, params){
#'   beta.hat <- coef(model)[2]
#'   beta.lwr <- confint(model)[2,1]
#'   beta.upr <- confint(model)[2,2]
#'   return( data.frame(beta=beta.hat, lwr=beta.lwr, upr=beta.upr))
#' }
#' 
#' Sims <- summarize.sims(Summary.Function, Params)
#' @export
summarize.sims <- function( summary.function, 
                            Params, sim.directory='./sim_directory', ...){

  target.dir <- paste(sim.directory, '/OutputFiles', sep='')
  files <- as.character( dir(path=target.dir ) )
	files <- data.frame(file = files,
	                    sim = get.sim.number(files),
	                    rep = get.rep.number(files))
	files <- files %>% arrange(sim, rep)
	
	out <- list()
	for(k in 1:nrow(files)){
	  i <- files$sim[k]
	  j <- files$rep[k]
		load(paste(target.dir,'/sim',i,'rep',j,'.RData', sep=''))
		temp1 <- summary.function(sim)
		if(is.vector(temp1)){
		  out[[k]] <- cbind( do.call(rbind, replicate(length(temp1), Params[i,], simplify=FALSE)), temp1 )
		}else if( is.matrix(temp1) | is.data.frame(temp1) ){
		  out[[k]] <- cbind( do.call(rbind, replicate(nrow(temp1),   Params[i,], simplify=FALSE)), temp1 )
		}
	}
	return( rbind_all(out) )
}



get.sim.number <- Vectorize(function( file ){
	temp <- strsplit(file, 'sim')
	temp <- strsplit(temp[[1]][2], 'rep')
	sim <- as.integer(temp[[1]][1])
	names(sim) <- NULL
	return(sim)
}, USE.NAMES=FALSE)

get.rep.number <- Vectorize(function( file ){
	temp <- strsplit(file, 'rep')
	temp <- strsplit(temp[[1]][2], '.RData')
	rep <- as.integer(temp[[1]][1])
	return(rep)
}, USE.NAMES=FALSE)
		

make.RNG.seeds <- function(seed, num.sims, num.reps){
  out <- list()
  set.seed(seed)
  for(i in 1:num.sims){
    for(j in 1:num.reps){
      out[[paste('sim',i,'rep',j,sep='')]] <- .Random.seed
      temp <- runif(10)
    }
  }
  return(out)
}

make.RNG.seed.names <- function(num.sims, num.reps){
  names <- expand.grid(sim=1:num.sims, rep=1:num.reps) %>%
    mutate(name = paste('sim',sim,'rep',rep,sep=''))
  return(names$name)
}
