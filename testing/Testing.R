library(dplyr)

user.function <- function(n, mu){
  return(data.frame(xbar = mean( rnorm(n,mean=mu))))
}

Params <- expand.grid(n  = c(100, 400),
                      mu = c(0, 10, 20))

function2 <- function(a,b){
  return(a^b)
}

x <- rnorm(5)

prep.sims(user.function, Params, 
          num.reps=40, sim.directory='~/Testing')

run.sims(sim.directory='~/Testing')

sims <- summarize.sims(identity, Params, '~/Testing')
sims %>% group_by(n,mu) %>%
  summarise(mean = mean(xbar),
            sd   = sd(xbar))
