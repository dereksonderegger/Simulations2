user.function <- function(n, mu){
  return(mean( rnorm(n,mean=mu) ) )
}

Params <- expand.grid(n  = c(100, 400),
                      mu = c(0, 10))

function2 <- function(a,b){
  return(a^b)
}

prep.sims(user.function, Params, 
          num.reps=10, sim.directory='~/Testing')

run.sims(sim.directory='~/Testing')
