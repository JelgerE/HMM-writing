library(gamlss)
library(gamlss.dist)
library(rgdal)
library(patchwork)

BEmu <- a / (a+b)
BEsigma <- (a*b) / ((a+b)^2 * (a+b+1))

BE(mu.link = 'logit', sigma.link = 'logit')

plot(dbeta(dat$SI[order(dat$SI)], shape1 = a, shape2 = b))


p <- ggplot(dat, aes(x = step))
p <- p + geom_line(aes(y = SI, colour = "state"))

# adding the relative humidity data, transformed to match roughly the range of the temperature
p <- p + geom_line(aes(y = rel_hum/5, colour = "Humidity"))

# now adding the secondary axis, following the example in the help file ?scale_y_continuous
# and, very important, reverting the above transformation
p <- p + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Relative humidity [%]"))

# modifying colours and theme options
p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Air temperature [°C]",
              x = "Date and time",
              colour = "Parameter")
p <- p + theme(legend.position = c(0.8, 0.9))
p


time <- ggplot(HMMdat2, aes(x=Time, y=state)) +
  geom_point()
time

shp <- readOGR(dsn = 'C:\\Users\\jelings\\OneDrive - UGent\\Documenten\\INBO\\Altusreid\\RStudio\\Altusried R\\shapefile\\new_river_shapefile.shp')
  
state1 <- ggplot(HMMdat2[HMMdat2$state == 1,] , aes(x=x, y=y)) + 
  ggtitle('State 1') +
  geom_point() +
  geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
  theme(plot.title=element_text(hjust=0.5))
state2 <- ggplot(HMMdat2[HMMdat2$state == 2,] , aes(x=x, y=y)) + 
  ggtitle('State 2') +
  geom_point() +
  geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
  theme(plot.title=element_text(hjust=0.5))
state1+state2
