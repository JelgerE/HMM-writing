library(moveHMM)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(depmixS4)
library(crawl)
library(EnvStats)
library(patchwork)
library(grDevices)

# Read data
fish <- read.csv(".\\data\\filtered_data_1_grayling_46909_HPE_RMSE_Vel.csv")
fish$Time <- as.POSIXct(fish$Time)

time_res <- c('5 secs', '10 secs', '15 secs', '30 secs', '60 secs', '2 mins', '5 mins')
fishreg <- crwMLE(data = fish, coord = c('UTM.x', 'UTM.y'), Time.name = c("Time") , time.scale = 'secs')

shp <- readOGR(dsn = 'C:\\Users\\jelings\\OneDrive - UGent\\Documenten\\INBO\\Altusreid\\RStudio\\Altusried R\\shapefile\\new_river_shapefile.shp')

for (res in time_res){
  # Construct movement model and regularize data
  fishregmod <- crwPredict(object.crwFit = fishreg, 
                           predTime = seq.POSIXt(from = min(fish$Time), to = max(fish$Time), by = res), 
                           return.type = 'minimal')
  
  # Isolate predicted positions
  fishpred <- fishregmod[fishregmod$locType == 'p', ]
  
  # Pred data for HMM(= calculate steps and TAs)
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
  #HMMmod <- depmix(step ~ 1, data = HMMdat2, nstates = 2,
  #                 family = Gamma(), respstart = c(2, 10))
  
  # Step length and TA
  #HMMmod <- depmix(list(step ~ 1,angle ~ 1), data = HMMdat2, nstates = 2,
  #                 family = list(gaussian(), gaussian()))
  
  # Step length and SI
  HMMmod <- depmix(list(step ~ 1,SI ~ 1), data = HMMdat2, 
                   family = list(Gamma(), gaussian()), respstart = c(5, 10, a, b, 1, 1), # ==> Why 6 values? 2 for steps, 2 for SI, 2 for????
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
              time.res = res,
              step.mean = mean(step), step.sd = sd(step),
              SI.mean = mean(SI), SI.sd = sd(SI),
              angle.mean = mean(angle), angle.sd = sd(angle))
  
  if (res == '5 secs') {
    mu_tot <- mu
  } else {
    mu_tot <- rbind(mu_tot, mu)
  }
  
  
  steps <- ggplot(HMMdat2, aes(x=step)) +
    geom_histogram(aes(y=..density..), #binwidth = (1/10),
                   colour = 1, fill = "white") +
    #geom_density() +
    geom_density(aes(y=..density.., color = state)) +
    #stat_function(fun = function(x) dgamma(x, shape = shape, scale = scale), linetype="dashed") + 
    #geom_vline(data=mu, aes(xintercept = step.mean, color=state)) + 
    ggtitle("Step Length") +
    theme(plot.title=element_text(hjust=0.5))
  
  SI <- ggplot(HMMdat2, aes(x=SI)) +
    geom_histogram(aes(y=..density..), #binwidth = (1/10),
                   colour = 1, fill = "white") +
    #geom_density() +
    geom_density(aes(y=..density.., color = state)) +
    #stat_function(fun = function(x) dgamma(x, shape = shape, scale = scale), linetype="dashed") + 
    #geom_vline(data=mu, aes(xintercept = step.mean, color=state)) + 
    ggtitle("Straightness Index") +
    theme(plot.title=element_text(hjust=0.5))
  
  TA <- ggplot(HMMdat2, aes(x=angle)) +
    geom_histogram(aes(y=..density..), #binwidth = (1/10),
                   colour = 1, fill = "white") +
    #geom_density() +
    geom_density(aes(y=..density.., color = state)) +
    #stat_function(fun = function(x) dgamma(x, shape = shape, scale = scale), linetype="dashed") + 
    #geom_vline(data=mu, aes(xintercept = step.mean, color=state)) + 
    ggtitle("Turning Angle") +
    theme(plot.title=element_text(hjust=0.5))
  
  jpeg(filename = paste('.\\figures\\', res, '_states.jpeg', sep=""))
  plot(steps + SI + TA +
         plot_layout(ncol = 2, guides = "collect") +
         plot_annotation(res, theme = theme(plot.title = element_text(hjust = 0.5))))
  dev.off()
  
  jpeg(filename = paste('.\\figures\\', res, '_map.jpeg', sep=""))
  if (res == "10 secs" | res == "30 secs" | res == "60 secs" | res == "5 mins"){
  state1 <- ggplot(HMMdat2[HMMdat2$state == 2,] , aes(x=x, y=y)) + 
    ggtitle("State 1") +
    geom_point() +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
    theme(plot.title=element_text(hjust=0.5))
  state2 <- ggplot(HMMdat2[HMMdat2$state == 1,] , aes(x=x, y=y)) + 
    ggtitle("State 2") +
    geom_point() +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
    theme(plot.title=element_text(hjust=0.5))
  }
  else {
    state1 <- ggplot(HMMdat2[HMMdat2$state == 1,] , aes(x=x, y=y)) + 
      ggtitle("State 1") +
      geom_point() +
      geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
      theme(plot.title=element_text(hjust=0.5))
    state2 <- ggplot(HMMdat2[HMMdat2$state == 2,] , aes(x=x, y=y)) + 
      ggtitle("State 2") +
      geom_point() +
      geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
      theme(plot.title=element_text(hjust=0.5))
  }
  plot(state1+state2+
         plot_annotation(res, theme = theme(plot.title = element_text(hjust = 0.5))))
  dev.off()
}

