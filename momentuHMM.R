library(crawl)
library(momentuHMM)
library(dplyr)

# Read data
fish <- read.csv(".\\data\\filtered_data_1_grayling_46909_HPE_RMSE_Vel.csv")
fish$Time <- as.POSIXct(fish$Time)

# Construct movement model and regularize data
fishreg <- crwMLE(data = fish, coord = c('UTM.x', 'UTM.y'), Time.name = c("Time") , time.scale = 'secs')

fishreg <- crwPredict(object.crwFit = fishreg, 
                      predTime = seq.POSIXt(from = min(fish$Time), to = max(fish$Time), by = '5 secs'), 
                      return.type = 'minimal')

# Isolate predicted positions
fishpred <- fishreg[fishreg$locType == 'p', ]

# Pred data for HMM(= calculate steps and TAs)
HMMdat <- moveHMM::prepData(fishpred[c('Time', 'mu.x', 'mu.y')], coordNames = c('mu.x', 'mu.y'), type = 'UTM')
#plot(HMMdat, compact=T)

# Calculate dStepLength and remove when dStep approaches 0 (indicating regularization in detection gaps)
# Breakpoint arbitrarily chosen based on rounding errors
HMMdat$dStep <- HMMdat$step - lag(HMMdat$step)
HMMdat <- HMMdat[(HMMdat$dStep > 3.172941e-08) | (HMMdat$dStep < -3.172941e-08), ]
HMMdat <- na.omit(HMMdat)
rownames(HMMdat) <- NULL 

# Calculate timesteps to identify where detections are removed (i.e. tracks break)
# Add in timestep for first value to run for-loop later
HMMdat$timestep <- HMMdat$Time - lag(HMMdat$Time)
HMMdat$timestep[1] <- as.difftime(5, format = '%S', units = 'secs')

# Construct trackno column and initialize with Track 1
HMMdat$track <- NA
HMMdat$track[1] <- 1

# Run for-loop to count tracks
for (n in 2:nrow(HMMdat)){
  if (HMMdat$timestep[n] > HMMdat$timestep[n-1]){
    HMMdat$track[n] <- HMMdat$track[n-1] + 1
  }
  else {
    HMMdat$track[n] <- HMMdat$track[n-1]
  }
}

# Get track lengths and remove tracks that are too short
Ltrack <- data.frame(table(HMMdat$track))
Ltrack <- Ltrack[Ltrack$Freq > 10,]
HMMdat <- HMMdat[HMMdat$track %in% Ltrack$Var1, ]

# Redefine trackno as animal ID for rerunning of prepData
HMMdat <- HMMdat[, c("x", "y", "Time", "track")]
colnames(HMMdat) <- c("x", "y", "Time", "ID")

# Re-run prepData for newly identified tracks
HMMdat2 <- momentuHMM::prepData(HMMdat, type='UTM')
#plot(HMMdat2, compact = T)

# Calculate Straightness Index (SI) over 30s
X <- split(HMMdat2, HMMdat2$ID)
for(d in 1:length(X)) {
  track <- data.frame(X[d])
  colnames(track) <- c("ID", "step", "angle", "Time", "x", "y")
  track$SI <- NA
  for (c in 5:(nrow(track)-1)){
    dX <- track$x[c] - track$x[c-4]
    dY <- track$y[c] - track$y[c-4]
    TDist <- sum(track$step[(c-4):c])
    LDist <- sqrt(dX^2 + dY^2)
    #print(LDist / TDist)
    track$SI[c] <- LDist / TDist
  }
  #print(colnames(track))
  if (d == 1) {
    HMMdat2 <- track
  }
  else {
    HMMdat2 <- rbind(HMMdat2, track)
  }
}

# Omit all NA values, reset index, and order by track
HMMdat2 <- na.omit(HMMdat2)
rownames(HMMdat2) <- NULL
HMMdat2 <- HMMdat2[order(as.numeric(rownames(HMMdat2))),,]

# Prep data for momentu package
momentuHMMdat <- momentuHMM::prepData(HMMdat2[, c('ID', 'Time', 'x', 'y', 'SI')], type='UTM')

# Set initial parameters for HMM (taken from depmix model)
stepPar0 <- c(0.16, 1.2, 0.13, 0.8)
anglePar0 <- c(0.1, 0.1)
dist <- list(stepDist='gamma', angleDist='wrpcauchy')

# Fit HMM
momentu_m <- momentuHMM::fitHMM(momentuHMMdat, nbStates=2, 
                                dist = list(step='gamma', angle='vm'),
                                Par0=list(step=stepPar0, angle=anglePar0))

# Print and plot HMM
momentu_m
#plot(momentu_m, compact=T)