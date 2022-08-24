# Draw histograms
# Run code in depmixS4.R first to get data

library(ggplot2)

ggplot(HMMdat2, aes(x=step)) +
  geom_histogram(aes(y=..density..), #binwidth = (1/10),
                 colour = 1, fill = 'white') +
  geom_density() +
  geom_density(aes(y=..density.., color = state)) +
  #stat_function(fun = function(x) dgamma(x, shape = shape, scale = scale), linetype='dashed') + 
  #geom_vline(data=mu, aes(xintercept = step.mean, color=state)) + 
  ggtitle('Histogram of step lengths') +
  theme(plot.title=element_text(hjust=0.5))

ggplot(HMMdat2, aes(x=SI)) +
  # Histogram of SI
  geom_histogram(aes(y=..density..), binwidth = (1/10),
                 color = 1, fill = 'white') +
  # Total density
  #geom_density(show.legend = T) +
  # Density curves per state
  geom_density(aes(y=..density.., color = state), show.legend = T) +
  # Beta distribution of SI
  #stat_function(fun = function(x) dbeta(x, shape1 = a, shape2 = b), linetype='dashed') + 
  #geom_vline(data=mu, aes(xintercept = SI.mean, color=state)) + 
  ggtitle('Histogram of SI')+
  theme(plot.title=element_text(hjust=0.5))

ggplot(HMMdat2, aes(x=angle, color = state)) +
  geom_histogram(aes(y=..density..), binwidth = (1/10),
                 colour = 1, fill = 'white') +
  geom_density()+
  #geom_vline(data=mu, aes(xintercept = angle.mean, color=state)) + 
  ggtitle('Histogram of angle')+
  theme(plot.title=element_text(hjust=0.5))
