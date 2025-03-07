# use "zika Brazil.rdata" from main data analysis

rm(list = ls())
library(SimInf)
library(Matrix)
library(truncnorm)
options(warn = 2)

#########
# Setup #
#########

# load data from https://github.com/jptrostle/SpatialSIRGPMC
load("zika Brazil data.rdata")

y <- cbind(zika_data[,, 1], zika_data[,, 2])[, 41:80]
rm(zika_data, states_names, area_km)
logdens <- log_density

# sim settings
ns <- 27
start_day <- 41
nt <- 40
start_inf <- 100
S0 <- 10
pop <- gridpop
# report rate is from main Zika analysis
report_rate <- c(0.01968, 0.00560, 0.00721, 0.02370, 0.01917, 0.00480, 0.00055, 0.00180, 0.00548, 0.00848, 0.04459, 0.00260, 0.00394, 0.01293, 0.00361, 0.00132, 0.00062, 0.00044, 0.01386, 0.00548, 0.00201, 0.00886, 0.00219, 0.00028, 0.00108, 0.00100, 0.00801)
this_gamma <- 1.2

nruns <- 10000 # should be around 100,000

##############################################################
# Functions From https://github.com/jptrostle/SpatialSIRGPMC #
##############################################################


# the Gillespie part
run_siminf <- function(ns, nt, start_day, grid_pop, s_gz, I0, logdens, thisbeta0, thisbeta1, thisphi, thisgamma, Adj) {
    
    nt <- nt + start_day
    
    thisbeta <- exp(thisbeta0 + thisbeta1 * logdens)
    thisbeta1 <- thisbeta[1]
    thisbeta2 <- thisbeta[2]
    thisbeta3 <- thisbeta[3]
    thisbeta4 <- thisbeta[4]
    thisbeta5 <- thisbeta[5]
    thisbeta6 <- thisbeta[6]
    thisbeta7 <- thisbeta[7]
    thisbeta8 <- thisbeta[8]
    thisbeta9 <- thisbeta[9]
    thisbeta10 <- thisbeta[10]
    thisbeta11 <- thisbeta[11]
    thisbeta12 <- thisbeta[12]
    thisbeta13 <- thisbeta[13]
    thisbeta14 <- thisbeta[14]
    thisbeta15 <- thisbeta[15]
    thisbeta16 <- thisbeta[16]
    thisbeta17 <- thisbeta[17]
    thisbeta18 <- thisbeta[18]
    thisbeta19 <- thisbeta[19]
    thisbeta20 <- thisbeta[20]
    thisbeta21 <- thisbeta[21]
    thisbeta22 <- thisbeta[22]
    thisbeta23 <- thisbeta[23]
    thisbeta24 <- thisbeta[24]
    thisbeta25 <- thisbeta[25]
    thisbeta26 <- thisbeta[26]
    thisbeta27 <- thisbeta[27]
    
    # list the adjacencies
    neighbors <- list()
    for (i in 1:ns) neighbors[[as.character(i)]] <- which(Adj[i, ] > 0)
    
    start_susc <- grid_pop
    start_inf <- rep(0, ns)
    start_inf[s_gz] <- I0
    start_susc <- start_susc - start_inf
    
    # local stuff initially
    transitions <- c()
    for (i in 1:ns) {
        thistextINF <- paste("S", i, " -> beta", i, "*S", i, "*I", i, "/(S", i, "+I", i, "+R", i, ") -> I", i, sep = "")
        thistextREC <- paste("I", i, " -> gamma*I", i, " -> R", i, sep = "")
        transitions <- c(transitions, thistextINF, thistextREC)
    }
    # spatial stuff
    for (i in 1:ns) {
        for (j in neighbors[[i]]) {
            # this is site i infects site j
            thistextINF <- paste("S", j, " -> phi*S", j, "*I", i, "/(S", j, "+I", j, "+R", j, ") -> I", j, sep = "")
            transitions <- c(transitions, thistextINF)
        }
    }
    
    compartments <- c()
    for (i in 1:ns) {
        compartments <- c(compartments, paste("S", i, sep = ""), paste("I", i, sep=""), paste("R", i, sep = ""))
    }
    
    u0 <- data.frame(S1 = rep(start_susc[1], 1), I1 = rep(start_inf[1], 1), R1 = rep(0, 1))
    for (i in 2:ns) {
        sCol <- paste("S", i, sep="")
        iCol <- paste("I", i, sep="")
        rCol <- paste("R", i, sep="")
        sAdd <- paste("u0$", sCol, " <- rep(start_susc[", i, "], 1)", sep = "")
        iAdd <- paste("u0$", iCol, " <- rep(start_inf[", i, "], 1)", sep = "")
        rAdd <- paste("u0$", rCol, " <- rep(0, 1)", sep="")
        eval(parse(text = sAdd))
        eval(parse(text = iAdd))
        eval(parse(text = rAdd))
    }
    
    model <- mparse(transitions = transitions,
                    compartments = compartments,
                    gdata = c(beta1 = thisbeta1,
                              beta2 = thisbeta2,
                              beta3 = thisbeta3,
                              beta4 = thisbeta4,
                              beta5 = thisbeta5,
                              beta6 = thisbeta6,
                              beta7 = thisbeta7,
                              beta8 = thisbeta8,
                              beta9 = thisbeta9,
                              beta10 = thisbeta10,
                              beta11 = thisbeta11,
                              beta12 = thisbeta12,
                              beta13 = thisbeta13,
                              beta14 = thisbeta14,
                              beta15 = thisbeta15,
                              beta16 = thisbeta16,
                              beta17 = thisbeta17,
                              beta18 = thisbeta18,
                              beta19 = thisbeta19,
                              beta20 = thisbeta20,
                              beta21 = thisbeta21,
                              beta22 = thisbeta22,
                              beta23 = thisbeta23,
                              beta24 = thisbeta24,
                              beta25 = thisbeta25,
                              beta26 = thisbeta26,
                              beta27 = thisbeta27,
                              phi = thisphi, gamma = thisgamma),
                    u0 = u0,
                    tspan = 1:nt)
    result <- run(model)
    simdata <- trajectory(result)

    return(simdata)
    
}

# the likelihood part
run_likelihood <- function(latent, ns, nt, start_day, p, nu) {
  
  dX <- matrix(0, nrow = ns, ncol = nt + start_day - 1)
  for (s in 1:ns) {
    this_col <- 3 * s
    dX[s, ] <- c(latent[1:(nt + start_day - 1), this_col] - latent[2:(nt + start_day), this_col])
  }
  
  dX <- c(dX) # to ensure properly vectorized
  y <- suppressWarnings(rnbinom(ns * (start_day + nt - 1), mu = p * dX, size = p * dX / (nu - 1)))
  y[is.na(y)] <- 0
  y <- matrix(y, nrow = ns, ncol = start_day + nt - 1, byrow = FALSE)
  return(y[, start_day:(start_day + nt - 1)])
  
}


##################################################################################
# Generating The Data from the Training Distribution to train the NN in VaNBayes #
##################################################################################

results = matrix(0, nrow = nruns, ncol = 1080+4) # 2187 from length of dataset, 3 from parameters

for (i in 1:nruns) {
    if (i %% 10 == 0) print(paste("On iter", i))
    
    # generate synthetic data set
    cand_beta0 <- runif(1, -3, 1)
    cand_beta1 <- runif(1, -1, 1)
    cand_phi <- rlnorm(1, -2, 1)
    cand_simdata <- run_siminf(ns, nt, start_day, pop, S0, start_inf, logdens, cand_beta0, cand_beta1, cand_phi, this_gamma, Adj)
    cand_nu <- runif(1, 1.01,10)
    cand_y <- run_likelihood(cand_simdata, ns, nt, start_day, report_rate, cand_nu)
    results[i,] = c(cand_beta0,cand_beta1,cand_phi,cand_nu,c(cand_y))
}

write.csv(results, file = "training_data.csv")
write.csv(c(y), file = "observed_data.csv")

# eof
