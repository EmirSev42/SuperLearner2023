# current location

### Breakdown.R ###
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Author: Emir S
## Created: Feb 6, 2023 
## Version: 
## Last-Updated: Jun 30 2023 (09:48) 
##           By:
##     Update #:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
### Commentary: 

# This is the informal "breakdown" of the main function (computing coefficients alpha);
# called "survSuperLearner" in the survSuperLearner package

# Here is Dr.Westling's Github with all of the codes:
# https://github.com/tedwestling/survSuperLearner/tree/master

# The "important" bit, where the coefficients are computed, is marked clealy as "THE IMPORTANT PART" with the comments
# the function that does this is called .survComputeCoef, and it is an internal function
# I hope this helps do a degree!
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# ----------------------------------------------------------------------- #
# STEP 1: THE INPUT WE WILL USE
# ----------------------------------------------------------------------- #

# this is identical to the example Westling uses on github
# we generate it simply so we can run the individual functions within the "survSuperLearner" function

# We have TWO predictors, X1 and X2, with time and event variables of appropriate dimension
# dimension of matrix X = [X1,X2]: n x 2
# dimension of time & event: n x 1

# We also set up the libraries to be used
# the library functions are individually very simple; they fit models for a given X
# and predict risks for a given newX for time (t)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
set.seed(92)
n <- 100
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))

S0 <- function(t, x) pexp(t, rate = exp(-2 + x[,1] - x[,2] + .5 * x[,1] * x[,2]), 
                          lower.tail = FALSE)
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))

G0 <- function(t, x) {
  as.numeric(t < 15) * .9 * pexp(t, rate = exp(-2 -.5 * x[,1] - .25 * x[,2] + .5 * x[,1] * x[,2]),
                                 lower.tail=FALSE)
}
C0 <- rbinom(n, 1, .1)
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C0 == 1] <- 0
C[C > 15] <- 15

time <- pmin(T, C)
event <- as.numeric(T <= C)

# each of these strings is the name of an individual function that computes risks for the respective method
# to see one for example:
# ?survSuperLearner::survSL.coxph()

event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.expreg", 
                                                "survSL.weibreg", "survSL.loglogreg", 
                                                "survSL.gam", "survSL.rfsrc")

# Let's also set Verbose=T, will see why in a second
verbose<-TRUE
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# ----------------------------------------------------------------------- #
# STEP 2: SETTING UP EVERYTHING THAT IS NEEDED FOR THE MAIN FUNCTION
# ----------------------------------------------------------------------- #

# THE FUNCTION CALL OF survSuperLearner LOOKS LIKE THIS:

# survSuperLearner(time,event,X,new.times,event.SL.library,verbose,cens.SL.library,newX,verbose)
# WHERE:

# time: the time variable
# event: the censoring indicator (1: dead, 0: alive)
# new.times: the times (t) to acquire risk predictions 
# event.SL.library: the libraries (base learners) to be used for the EVENT calculation
# cens.SL.library: the libraries (base learners) to be used for the CENSORING calculation
# newX: the "test data" or matrix to predict the risks for
# verbose: details or not

# NOTE: it seems that generally event.SL.library = cens.SL.library, though the function does not demand this

# So now we see what the function does with the inputs provided, and hopefully we will get to how the coefficients are calculated intenrally

# first save time and event as numbers
time <- as.numeric(time)
event <- as.numeric(event)

# weights if specified, we set to null
obsWeights<-NULL
# newX if specified, we set to null
newX<-NULL

# set observation weights to 1 if not provided
if (is.null(obsWeights)) obsWeights <- rep(1, length(time))
if (is.null(newX)) newX <- X

# set up variralble names, with N and p equal to the row and
# column numbers of the matrix X
varNames <- colnames(X)
N <- dim(X)[1L]
p <- dim(X)[2L]

# controls set to empty lists
control<-list()
cvControl<-list()

# CREATE LIBRARY FUNCTION
# this function creates the libraries to be used as a list object to be referenced later

.createLibrary <- function (survSL.library)  {
  if (is.character(survSL.library)) {
    k <- length(survSL.library)
    whichScreen <- matrix(1, nrow = 1, ncol = k)
    screenAlgorithm <- "All"
    library <- data.frame(predAlgorithm = survSL.library, rowScreen = 1,
                          stringsAsFactors = FALSE)
  }
  else if (is.list(survSL.library)) {
    predNames <- sapply(survSL.library, FUN = "[", 1)
    NumberScreen <- (sapply(survSL.library, FUN = length) - 1)
    if (sum(NumberScreen == 0) > 0) {
      for (ii in which(NumberScreen == 0)) {
        SL.library[[ii]] <- c(SL.library[[ii]], "All")
        NumberScreen[ii] <- 1
      }
    }
    screenAlgorithmFull <- unlist(lapply(survSL.library, FUN = "[", -1))
    screenAlgorithm <- unique(screenAlgorithmFull)
    library <- data.frame(predAlgorithm = rep(predNames,
                                              times = NumberScreen), rowScreen = match(screenAlgorithmFull,
                                                                                       screenAlgorithm), stringsAsFactors = FALSE)
  }
  else {
    stop("format for survSL.library is not recognized")
  }
  out <- list(library = library, screenAlgorithm = screenAlgorithm)
  return(out)
}

if (is.null(control$event.t.grid)) {
  control$event.t.grid <- seq(0, max(time[event == 1]), length.out = 250)
} else {
  if (!is.null(dim(control$event.t.grid)) && ncol(control$event.t.grid) > 0) stop("t.grid must be an (k x 1)  numeric vector.")
  control$event.t.grid <- sort(unique(as.numeric(control$event.t.grid)))
  if(any(is.na(control$event.t.grid))) stop("No missing values allowed in t.grid")
  if (any(control$event.t.grid < 0)) stop("Values in t.grid must be non-negative")
  if (any(control$event.t.grid > max(time))) stop("Values in t.grid must be no large than max(time)")
}

if (is.null(control$cens.t.grid)) {
  control$cens.t.grid <- seq(0, max(time[event == 0]), length.out = 250)
  control$cens.t.grid <- control$cens.t.grid - max(min(diff(sort(unique(time)))) / 2, 1e-5)
  control$cens.t.grid <- c(0, control$cens.t.grid[control$cens.t.grid > 0])
} else {
  if (!is.null(dim(control$cens.t.grid)) && ncol(control$cens.t.grid) > 0) stop("t.grid must be an (k x 1)  numeric vector.")
  control$cens.t.grid <- sort(unique(as.numeric(control$cens.t.grid)))
  if(any(is.na(control$cens.t.grid))) stop("No missing values allowed in t.grid")
  if (any(control$cens.t.grid < 0)) stop("Values in t.grid must be non-negative")
  if (any(control$cens.t.grid > max(time))) stop("Values in t.grid must be no large than max(time)")
}

# these cretae control parameters for the algorithm
# i.e:


# (THIS IS THE TIME GRID WE DISCUSSED IN OUR MEETING!) --- ------ ------ ---

# event.t.grid: Grid of times to use to approximate the integral in the risk function for the conditional survival function of the event. 
# Defaults to 250 points equally spaced between 0 and the last uncensored follow-up time.

# initWeightAlg: initial weighting algorithm defaulting to "survSL.rfsrc"; RANDOM FOREST

# for more details:
# ?survSuperLearner.control
# ?survSuperLearner.CV.control
control <- do.call("survSuperLearner.control", control)
cvControl <- do.call("survSuperLearner.CV.control", cvControl)

# call the create library function
event.library <- .createLibrary(event.SL.library)
cens.library <- .createLibrary(cens.SL.library)

# additional setups and checks
event.k <- nrow(event.library$library)
cens.k <- nrow(cens.library$library)
event.kScreen <- length(event.library$screenAlgorithm)
cens.kScreen <- length(cens.library$screenAlgorithm)
event.Z <- array(NA, dim = c(N, length(control$event.t.grid), event.k))
cens.Z <- array(NA, dim = c(N, length(control$cens.t.grid), cens.k))
event.libraryNames <- cbind(predAlgorithm = event.library$library$predAlgorithm,
                            screenAlgorithm = event.library$screenAlgorithm[event.library$library$rowScreen])
cens.libraryNames <- cbind(predAlgorithm = cens.library$library$predAlgorithm,
                           screenAlgorithm = cens.library$screenAlgorithm[cens.library$library$rowScreen])
if (p < 2 & !identical(event.library$screenAlgorithm, "All")) {
  warning("Screening algorithms specified in combination with single-column X.")
}

if (p < 2 & !identical(cens.library$screenAlgorithm, "All")) {
  warning("Screening algorithms specified in combination with single-column X.")
}

# catch if libraries give error
event.errorsInCVLibrary <- event.errorsInLibrary <- rep(0, event.k)
cens.errorsInCVLibrary <- cens.errorsInLibrary <- rep(0, cens.k)

# cv folds function
# this creates folds to be used for the cross-validation calculated library predictions

.survCVFolds <- function (N, id, event, cvControl) {
  if (!is.null(cvControl$validRows)) return(cvControl$validRows)
  stratifyCV <- cvControl$stratifyCV
  shuffle <- cvControl$shuffle
  V <- cvControl$V
  if (!stratifyCV) {
    if (shuffle) {
      if (is.null(id)) {
        validRows <- split(sample(1:N), rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(sample(1:n.id), rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
    else {
      if (is.null(id)) {
        validRows <- split(1:N, rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(1:n.id, rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
  }
  else {
    # if (sum(event) < V | sum(1-event) < V) {
    #   stop("number of (event = 1) or (event = 0) is less than the number of folds")
    # }
    if (shuffle) {
      if (is.null(id)) {
        event.0 <- which(event == 0)
        event.1 <- which(event == 1)
        rows.0 <- split(sample(event.0), rep(1:V, length = length(event.0)))
        rows.1 <- split(sample(event.1), rep(1:V, length = length(event.1)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          if (length(rows.0) >= vv) {
            if (length(rows.1) >= vv) validRows[[vv]] <- c(rows.0[[vv]], rows.1[[vv]])
            else validRows[[vv]] <- rows.0[[vv]]
          } else {
            validRows[[vv]] <- rows.1[[vv]]
          }
        }
      }
      else {
        stop("Stratified sampling with id not currently implemented. Either remove id or set control(stratifyCV = FALSE).")
      }
    }
    else {
      if (is.null(id)) {
        within.split <- suppressWarnings(tapply(1:N,
                                                INDEX = event, FUN = split, rep(1:V)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          validRows[[vv]] <- c(within.split[[1]][[vv]],
                               within.split[[2]][[vv]])
        }
      }
      else {
        stop("Stratified sampling with id not currently implemented. Either remove id or set control(stratifyCV = FALSE).")
      }
    }
  }
  return(validRows)
}

# call above function
validRows <- .survCVFolds(N = N, id = NULL, event = event, cvControl = cvControl)
if (is.null(id)) id <- seq(N)
id <- seq(N)

# function to get cross-validated survivals
.crossValFUN <- function(valid, time, event, dataX, id, obsWeights, t.grid,
                         library, kScreen, k, p) {
  foldNum <- as.numeric(which(unlist(lapply(validRows, function(v) all.equal(v, valid) == TRUE))))
  tempLearn <- dataX[-valid, , drop = FALSE]
  tempTime <- time[-valid]
  tempEvent <- event[-valid]
  tempValid <- dataX[valid, , drop = FALSE]
  tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
  tempId <- id[-valid]
  tempObsWeights <- obsWeights[-valid]
  for (s in seq(kScreen)) {
    if(verbose) message(paste("CV ", library$screenAlgorithm[s],
                              ", fold ", foldNum, sep = ""))
    screen_fn <- get(library$screenAlgorithm[s])
    testScreen <- try(do.call(screen_fn, list(time = tempTime,
                                              event = tempEvent,
                                              X = tempLearn, id = tempId,
                                              obsWeights = tempObsWeights)))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,",
                    library$screenAlgorithm[s], ", with All()",
                    "\n "))
      tempWhichScreen[s, ] <- TRUE
    }
    else {
      tempWhichScreen[s, ] <- testScreen
    }
    if (verbose) {
      message(paste("Number of covariates in ", library$screenAlgorithm[s],
                    " is: ", sum(tempWhichScreen[s, ]), sep = ""))
    }
  }
  
  uniqueScreen <- unique(tempWhichScreen)
  screenMap <- apply(uniqueScreen, 1, function(row) which(apply(tempWhichScreen, 1, function(row2) all.equal(row, row2) == TRUE)))
  
  out <- array(NA, dim = c(nrow(tempValid), length(t.grid), k))
  
  for (predAlg in unique(library$library$predAlgorithm)) {
    if (verbose) message(paste("CV ", predAlg, ", fold ", foldNum, sep = ""))
    pred_fn <- get(predAlg)
    for(j in seq(nrow(uniqueScreen))) {
      testAlg <- try(do.call(pred_fn, list(time = tempTime, event = tempEvent,
                                           X = subset(tempLearn, select = uniqueScreen[j,], drop = FALSE),
                                           newX = subset(tempValid, select = uniqueScreen[j,], drop = FALSE),
                                           new.times = t.grid,
                                           id = tempId,
                                           obsWeights = tempObsWeights)))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", predAlg,
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      } else {
        libraryRows <- which(library$library$predAlgorithm == predAlg & library$library$rowScreen %in% unlist(screenMap[j]))
        for (row in libraryRows) {
          out[,,row] <- testAlg$pred
        }
      }
    }
  }
  
  
  invisible(list(out = out))
}

# Here we compute risk predictions for each library in each of the folds
# So we have cross-validated risk and cens predictions for each library
                                                                
#  apply the cross validation function to the folds, set the event = event (for event probability calculations)                                     
time_train_start <- proc.time()
event.crossValFUN_out <- lapply(validRows, FUN = .crossValFUN,
                                time = time, event = event, dataX = X, id = id, obsWeights = obsWeights,
                                t.grid = control$event.t.grid, library = event.library,
                                kScreen = event.kScreen, k = event.k, p = p)


# catch errors here
for(v in 1:cvControl$V) {
  event.Z[validRows[[v]], ,] <- event.crossValFUN_out[[v]]$out
}

event.errorsInCVLibrary <- apply(event.Z, 3, anyNA)
if (sum(event.errorsInCVLibrary) > 0) event.Z[, , as.logical(event.errorsInCVLibrary)] <- 0
if (all(event.Z == 0))  stop("All algorithms dropped from event library")

# same function call as earlier, but this time with event = 1-event, so we're calculating censoring probabilities instead of risk
cens.crossValFUN_out <- lapply(validRows, FUN = .crossValFUN,
                               time = time, event = 1-event, dataX = X, id = id, obsWeights = obsWeights,
                               t.grid = control$cens.t.grid, library = cens.library,
                               kScreen = cens.kScreen, k = cens.k, p = p)

# catch errors again
for(v in 1:cvControl$V) {
  cens.Z[validRows[[v]], ,] <- cens.crossValFUN_out[[v]]$out
}

cens.errorsInCVLibrary <- apply(cens.Z, 3, anyNA)
if (sum(cens.errorsInCVLibrary) > 0) cens.Z[, , as.logical(cens.errorsInCVLibrary)] <- 0
if (all(cens.Z == 0))  stop("All algorithms dropped from censoring library")



# ----------------------------------------------------------------------- #
# STEP 3: THE IMPORTANT PART
# ----------------------------------------------------------------------- #

# Here we try to understand a single step of the algorithm used to compute the weights, alpha
# it is done via this function, using it iteratively
# before we iterate, let's see what it does for one step                                                                

.survcomputeCoef <- function(time, event, t.vals, cens.vals, preds, obsWeights) {
  if(ncol(preds) == 1) return(1)
  cens.vals[cens.vals < 1e-4] <- 1e-4
  out <- 1 - as.numeric(time <= t.vals) * event / cens.vals
  fit.nnls <- nnls::nnls(sqrt(obsWeights) * preds, sqrt(obsWeights) * out)
  # ind <- as.numeric(time > t.vals)
  # obsweight <- obsWeights * event / cens.vals
  # fit.nnls <- nnls::nnls(sqrt(obsweight) * preds, sqrt(obsweight) * ind)
  coef <- coef(fit.nnls)
  if(sum(coef) == 0) {
    warning("All coefficients in NNLS fit are zero.")
    coef <- rep(1,length(coef))
  }
  coef  / sum(coef)
}

# NOTE: the function inputs are deceptive! It's not the "time" and "event" that we set up in Part 1
# We will set this up now.

# the k's: event dimensioms
event.k <- dim(event.Z)[3]
cens.k <- dim(cens.Z)[3]

# N: length of the actual data (the time value we supplied to the function)
N <- length(time)
                                                                
# length of time for event and cens: the length of the grid. Length is always 250
event.n.time <- length(control$event.t.grid)
cens.n.time <- length(control$cens.t.grid)

# epsilon: lagged differences of unique time values are calculated and it's minimum is taken
# that or 1e-5, whichever is larger. I'm not sure how this is justified.
epsilon <- max(min(diff(sort(unique(time)))), 1e-5)

# event.Z.long and cens.Z.long: these are predictions from base learners in long form, survival and censoring prob. respectively
# The dimensions are (N*250 x no. libraries)
# i.e one prediction for each time grid point, from t1 to t250, repeated N times, per column (per library)
                                                                
event.Z.long <- matrix(NA, nrow = N * event.n.time, ncol = event.k)
for(j in seq(dim(event.Z)[3])) event.Z.long[,j] <- c(event.Z[,,j])
cens.Z.long <- matrix(NA, nrow = N * cens.n.time, ncol = cens.k)
for(j in seq(dim(cens.Z)[3])) cens.Z.long[,j] <- c(cens.Z[,,j])

# event.Z.obs and cens.Z.obs: basically the as above BUT instead of a prediction for each time grid point,
# the predictions are for the times ACTUALLY OBSERVED; i.e the time input to the main function
# these are computed via function interpolation; where the Y is the long form survival probs above, and
# X is the inputted time

event.Z.obs <- matrix(NA, nrow = N, ncol = event.k)
for(i in seq(N)) {
  for(j in seq(event.k)) {
    event.Z.obs[i,j] <- stats::approx(control$event.t.grid, event.Z[i,,j], xout = time[i], method = 'constant', rule = 2)$y
  }
}

cens.Z.obs <- matrix(NA, nrow = N, ncol = cens.k)
for(i in seq(N)) {
  for(j in seq(cens.k)) {
    cens.Z.obs[,j] <- stats::approx(c(-1,control$cens.t.grid), c(1,cens.Z[i,,j]), xout = time[i] - epsilon, method = 'constant', rule = 2)$y
  }
}


# next few are simple; the weights, observed time, observed event and the time grid, all in long form (repeated N times)
                                                                
obsWeights.event.long <- rep(obsWeights, event.n.time)
obsWeights.cens.long <- rep(obsWeights, cens.n.time)
time.event.long <- rep(time, event.n.time)
time.cens.long <- rep(time, cens.n.time)
event.event.long <- rep(event, event.n.time)
event.cens.long <- rep(event, cens.n.time)
event.t.grid.long <- rep(control$event.t.grid, each = N)
cens.t.grid.long <- rep(control$cens.t.grid, each = N)

# next: an initial weighting algorithm is set.
# in this case, it is just a random forest, set on the CENSORING as the event
# that is: event = 1 - event                                                              
# the new time is our time input minus epsilon
initWeightAlg <- get(control$initWeightAlg)
initFit <- initWeightAlg(time = time, event = 1 - event, X = X, newX = X,
                         new.times = time - epsilon,
                         obsWeights = obsWeights, id = id)
                                                                
# the DIAGONAL of the predictions are taken
# so each prediction is for the corresponding observed time.
# i.e: for i=1 take t1-eps, i=2 take e2-eps, i=3 take t3-eps etc...
# then put these into long form (replicate)
obs.cens.vals <- rep(diag(initFit$pred), length(control$event.t.grid))

# next, we get an initial estimate of coefficients
# set the coefficients as zero (length is equal to the number of candidate learners)
S.coef <- rep(0, event.k)

# ----------- WE NOW HAVE ALL THE INPUTS WE NEED FOR THE .suvComputeCoef function! ------------ #
# But how does it work?

 .survcomputeCoef <- function(time, event, t.vals, cens.vals, preds, obsWeights) {
   # IF only one set of predictions exist (so there is only one library), then we return 1 for all coefs
  if(ncol(preds) == 1) return(1)
   # next: for any censoring values too small, set them to 1e-4
  cens.vals[cens.vals < 1e-4] <- 1e-4
   # next: set out = 1 - { I(time < TimeGrid) * Event / CensoringValues }
   # where ALL ARE IN LONG FORM AS SET UP EARLIER
  out <- 1 - as.numeric(time <= t.vals) * event / cens.vals
  # next: use NNLS (non-negative least squares) to solve min ||Ax-b||_2
  # WHERE A = sqrt(obsWeights) * preds, and B = sqrt(obsWeights) * out
  fit.nnls <- nnls::nnls(sqrt(obsWeights) * preds, sqrt(obsWeights) * out)
  # grab the coefficients of this optimization
  coef <- coef(fit.nnls)
   # warn if all are zero
  if(sum(coef) == 0) {
    warning("All coefficients in NNLS fit are zero.")
    coef <- rep(1,length(coef))
  }
   # take them as proportions; i.e standardize
  coef  / sum(coef)
}                                                               

# first this function is applied once to grab the initial coefficients for the SURVIVAL curve (S)
# # what's important for now are the INPUTS:
# time, event, tvals,censvals, obsweights and preds are all LONG FORM                                                          
S.coef[!event.errorsInLibrary] <- .survcomputeCoef(time = time.event.long,
                                                   event = event.event.long,
                                                   t.vals = event.t.grid.long,
                                                   cens.vals = obs.cens.vals,
                                                   preds = event.Z.long[,!event.errorsInLibrary, drop=FALSE],
                                                   obsWeights = obsWeights.event.long)

# next, this is created
obs.event.vals <- rep(c(event.Z.obs %*% S.coef), length(control$cens.t.grid))                                                                
                                                                
# here is what this is:
# take the predictions for the observed values from the base learners:
event.Z.obs%>%head()
event.Z.obs%>%dim() 

 # multiply it by the initial coefficient estimate
c(event.Z.obs %*% S.coef)%>%head()
c(event.Z.obs %*% S.coef)%>%length()

 # then long form it (replicate tgrid = 250 times). Here is how they look:
obs.event.vals%>%head()
obs.event.vals%>%length()                                                               

 # Now, how does the algorithm loop?
# ----------------------------------------------------------------------- #
# STEP 4: ONE STEP OF THE LOOP
# ----------------------------------------------------------------------- #
                                                                
# first, if obs.cens.vals and obs.event.vals are not empty, set the old to be the current versions
if(!is.null(obs.cens.vals)) obs.cens.vals.old <- obs.cens.vals
if(!is.null(obs.event.vals)) obs.event.vals.old <- obs.event.vals

# next, update the obs.cens.vals by using the compute coef function on the CENSORING
G.coef <- rep(0, cens.k)
G.coef[!cens.errorsInLibrary] <- .survcomputeCoef(time = time.cens.long,
                                                  event = 1 - event.cens.long,
                                                  t.vals = cens.t.grid.long,
                                                  cens.vals = obs.event.vals,
                                                  preds = cens.Z.long[,!cens.errorsInLibrary, drop=FALSE],
                                                  obsWeights = obsWeights.cens.long)

# multiply observed censoring probabilities by the computed coefficient for G
# and then long form it
obs.cens.vals <- rep(c(cens.Z.obs %*% G.coef), length(control$event.t.grid))
# if there are any zeroes, set them to the minimum non-zero value
obs.cens.vals[obs.cens.vals == 0] <- min(obs.cens.vals[obs.cens.vals > 0])

# old values
obs.cens.vals.old%>%head()
# new values
obs.cens.vals%>%head()

# update observed event vals similarly, using the EVENT as the input
S.coef[!event.errorsInLibrary] <- .survcomputeCoef(time = time.event.long, 
                                                   event = event.event.long,
                                                   t.vals = event.t.grid.long, 
                                                   cens.vals = obs.cens.vals,
                                                   preds = event.Z.long[,!event.errorsInLibrary, drop=FALSE],
                                                   obsWeights = obsWeights.event.long)

# update the event values by multiplying the predictions for observed times
# by the current coefficients for S, and log form it
obs.event.vals <- rep(c(event.Z.obs %*% S.coef), length(control$cens.t.grid))
# similarly replace zeroes
obs.event.vals[obs.event.vals == 0] <- min(obs.event.vals[obs.event.vals > 0])

# set the differences between current and old values
cens.delta <- max(abs(obs.cens.vals - obs.cens.vals.old))
cens.delta
event.delta <- max(abs(obs.event.vals - obs.event.vals.old))
event.delta

# check convergence condition
cens.delta + event.delta < 1e-5

# if it didnt converge, go to the next iteration; so run .survComputeCoef for both the event and censoring values
# using the current "obs.cens.vals" value, and grab the new coefficients...                                                                

                                                                
