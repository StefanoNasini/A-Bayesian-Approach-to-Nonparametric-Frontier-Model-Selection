

rm(list=ls(all=TRUE))
gc(reset = TRUE)

library(lpSolve)
library('data.table')
library(truncnorm)

# parallel ----------------------------------------------------------------
# %dopar%
library(doFuture)
library(foreach)
library(doRNG)
library(progressr)
library(progress)

myCores <- 2 # more than this does not work for large number of replications...
myParallel <- TRUE



###############################################
###############################################

source("Bayesian_DEA_MCMC_functions.R")


###############################################
###############################################


#---------------------------------------------
# Generate some data
#---------------------------------------------

data <- data.table(read.table("mkt2015-data.txt", header=TRUE))
head(data)

myYear <- 2001
#    Year 2001 ---------------------------------------------------------------

XX <- as.matrix(data[data$year == myYear, paste0("x",1:5)])
YY <- as.matrix(data[data$year == myYear, paste0("y",1:5)])
WW <- as.matrix(data[data$year == myYear, paste0("w",1:5)])

head(XX)

round( cor(XX), 2)

round( cor(YY), 2)
round( cor(WW), 2)


k = nrow(XX)

P_I = WW

m = ncol(YY)
n = ncol(XX)

Proposal_var = 0.1

getwd()

# change to the directory with the 'year'
setwd(paste0(getwd(),"/CRS/",myYear))

pbar <-
    progress_bar$new(
        format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
        total = k,
        complete = "=",   # Completion bar character
        incomplete = "-", # Incomplete bar character
        current = ">",    # Current bar character
        clear = FALSE,    # If TRUE, clears the bar when finish
        width = 100,      # Width of the progress bar
        force = TRUE
    )

cat(paste0("Year is ",myYear))

### setup parallel

parallel <- myParallel
cores <- myCores

handlers('progress')
handlers(global = TRUE)
handlers(
    list(
        handler_progress(
            format   = ":spin :message [:bar] Done :current/:total (:percent) in :elapsed ETA: :eta",
            width    = 60,
            complete = "#"
        )
    ),
    append = TRUE
)

registerDoFuture()
if (parallel == TRUE & cores == 1) {
    parallel <- FALSE
}
if (parallel == TRUE & cores > 1) {
    plan(multisession, workers = cores)
} else {
    plan(sequential)
}

runRange <- 1:k

# for(kkk in 1:k)
progressr::with_progress(enable = TRUE, {
    pprogress <- progressr::progressor(along = runRange)
    foreach(kkk = runRange) %dopar%
        {
            out = MetropolisHansongs_z_posterior_CRS(kk = kkk, tolerance= 1.0e-07, ITER_max= 30000, ITER_min= 200, proposal_var = Proposal_var, technology = 0, H = 1)
            thisDate <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
            pprogress(paste0("(",thisDate,")"))
        }
}, handlers = progressr::handlers("progress"))





