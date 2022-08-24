library(moveHMM)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(depmixS4)
library(crawl)
library(EnvStats)
library(sf)

# Read data
fish <- read.csv(".\\data\\Grayling_46909.csv")
fish$dt <- as.POSIXct(fish$dt)
names(fish)[names(fish) == 'dt'] <- 'Time'
fish <- fish %>% select(Time, lon, lat)

fish <- st_as_sf(fish, coords = c('lon', 'lat'), crs='epsg:4326')
fish <- st_transform(fish, crs = 'epsg:32632')

# Construct movement model and regularize data
fishreg <- crwMLE(data = fish, Time.name = c("Time") , time.scale = 'secs')

fishreg <- crwPredict(object.crwFit = fishreg, 
                      predTime = seq.POSIXt(from = min(fish$Time), to = max(fish$Time), by = '5 secs'), 
                      return.type = 'minimal')

# Isolate predicted positions
fishpred <- fishreg[fishreg$locType == 'p', ]

# Prep data for HMM(= calculate steps and TAs)
HMMdat <- moveHMM::prepData(fishpred[c('Time', 'mu.x', 'mu.y')], coordNames = c('mu.x', 'mu.y'), type = 'UTM')
#plot(HMMdat, compact=T)

# Calculate dStepLength and remove when dStep approaches 0 (indicating regularization in detection gaps)
# Breakpoint arbitrarily chosen based on rounding errors
HMMdat$dStep <- HMMdat$step - lag(HMMdat$step)
HMMdat <- HMMdat[(HMMdat$dStep > 3.172941e-08) | (HMMdat$dStep < -3.172941e-08), ]
HMMdat <- drop_na(HMMdat)
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
HMMdat2 <- moveHMM::prepData(HMMdat, type='UTM')
#plot(HMMdat2, compact = T)

# Calculate Straightness Index (SI) over 30s
X <- split(HMMdat2, HMMdat2$ID)
for(d in 1:length(X)) {
  track <- data.frame(X[d])
  colnames(track) <- c("ID", "step", "angle", "x", "y", "Time")
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

HMMdat2 <- na.omit(HMMdat2)
rownames(HMMdat2) <- NULL
HMMdat2 <- HMMdat2[order(as.numeric(rownames(HMMdat2))),,]

# Get tracklengths to use in HMM later
tracklengths <- data.frame(table(as.numeric(HMMdat2$ID)))

# Find data distributions
# Step length
gamma <- egamma(HMMdat2$step)
shape <- as.numeric(gamma$parameters[1])
scale <- as.numeric(gamma$parameters[2])

# Straightness Index
beta <- ebeta(HMMdat2$SI)
a <- as.numeric(beta$parameters[1])
b <- as.numeric(beta$parameters[2])

# Construct HMMs using depmixS4
# Just step length
#HMMmod <- depmix(step ~ 1, data = HMMdat2, nstates = 2)

# Step length and TA
#HMMmod <- depmix(list(step ~ 1,angle ~ 1), data = HMMdat2, nstates = 2,
#                 family = list(gaussian(), gaussian()))

# Step length and SI
HMMmod <- depmix(list(step ~ 1,SI ~ 1), data = HMMdat2, 
                 family = list(Gamma(), gaussian()), respstart = c(5, 10, a, b, 1, 1),
                 nstates = 2, ntimes = tracklengths$Freq)


# Fit HMM to data and calculate states+probabilities
fitHMMmod <- fit(HMMmod) 
# ==> Warning message: In em.depmix(object = object, maxit = emcontrol$maxit, tol = emcontrol$tol,  :
                       #Log likelihood decreased on iteration 33 from -6417.22339094861 to -6417.61990695331
HMMdat2$state <- as.factor(depmixS4::viterbi(fitHMMmod)$state)
HMMdat2$prob1 <- depmixS4::viterbi(fitHMMmod)$S1
HMMdat2$prob2 <- depmixS4::viterbi(fitHMMmod)$S2

# Get mean values for states
mu <- ddply(HMMdat2, 'state', summarise, 
            step.mean = mean(step), step.sd = sd(step),
            SI.mean = mean(SI), SI.sd = sd(SI),
            angle.mean = mean(angle), angle.sd = sd(angle))


