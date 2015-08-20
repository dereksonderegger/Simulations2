#' @export
prep.sims <- function( sim.function, param.matrix, 
                       num.reps=10, sim.directory='sim_directory'){

	num.sims <- dim(param.matrix)[1]

	# create the sim_directory if necessary
	if( !file.exists(sim.directory) ){ 
 		dir.create(sim.directory)
	  dir.create(paste(sim.directory,'/InputFiles', sep=''))
	  dir.create(paste(sim.directory,'/OutputFiles',sep=''))
	}

	# I should figure out how to grab the global envirnment, add 
	# a link to the sim.function and param.matrix and save the 
	# whole she-bang
	Sim.Env <- globalenv()
	Sim.Env$Sim.Function <- sim.function
	for( i in 1:num.sims ){
	  Sim.Env$Params  <- param.matrix[i,]
	  for(j in 1:num.reps){
	    Sim.Env$Output.File <- paste(sim.directory,'/OutputFiles/sim',i,'rep',j,'.RData',sep='')
	    save(list  = ls(all.names=TRUE, envir=Sim.Env),
	         file  = paste(sim.directory,'/InputFiles','/sim',i,'rep',j,'.RData',sep=''),
	         envir = Sim.Env )
	  }
	}
	
	# Generate a file that contains one command for each file
	file.create(paste(sim.directory,'/CommandFile.txt',sep=''))
	CommandFile <- file(paste(sim.directory,'/CommandFile.txt',sep=''), open='a' )
	script <- system.file('AnalyzeOneFile.R', package='Simulations2')
	for(i in 1:num.sims){
	  for(j in 1:num.reps){
	    input.file <- paste(sim.directory,'/InputFiles','/sim',i,'rep',j,'.RData',sep='')
	    writeLines(paste('Rscript ', script, ' ', input.file), CommandFile)
	  }
	}
	close(CommandFile)
} 

#' @export
run.sims <- function(sim.directory){
  files <- list.files(paste(sim.directory,'/InputFiles',sep=''))
  script <- system.file('AnalyzeOneFile.R', package='Simulations2')
  foreach(file = files) %dopar% {
    foo <- system( paste('Rscript ', script, ' ', sim.directory,'/InputFiles/',file, sep='') ) 
  }
  return(invisible(TRUE))
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

  files <- as.character( dir(path=sim.directory) )
	files <- data.frame(file = files,
	                    sim = get.sim.number(files),
	                    rep = get.rep.number(files))
	files <- files %>% arrange(sim, rep)
	
	out <- list()
	for(k in 1:nrow(files)){
	  i <- files$sim[k]
	  j <- files$rep[k]
		load(paste(sim.directory,'/sim',i,'rep',j,'.RData', sep=''))
		temp1 <- summary.function(sim)
		temp2 <- cbind( do.call(rbind, replicate(nrow(temp1), Params[i,], simplify = FALSE)),
		                 temp1 )
		out[[k]] <- temp2  
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
		
