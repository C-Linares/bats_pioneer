# Script: Ploting_pop_models.R
# The purpose of this script is to replicate what Dr.Cruz did in here script :
# https://github.com/quantitativeconservationlab/AppPopnEco/blob/master/CountBayesPlotting.R
# 
# 
# Carlos Linares Dec 2023.
# 
# load libraries.
# load packages:
library( tidyverse ) 
options( dplyr.width = Inf, dplyr.print_min = 100 ) #all columns are printed and more rows are displayed
library( jagsUI ) 


# Load data 
# no data here apparently

# get summary of the model from pop_models.R

summary(m1)

# define result to plot

mr<-m1 # why doingthi?

############## trace plots ############
#plot( mr ) #plots traces and posterior densities for all parameters
par( mfrow = c( 2, 2 ), ask = F, mar = c(3,4,2,2) )
#detection parameters
#intercept 
traceplot( mr, parameters = c( 'int.det') )
#random intercepts for observer
traceplot( mr, parameters = c( 'eps.det') )
#fixed effect
traceplot( mr, parameters = c( 'alpha') )
#for abundance
traceplot( mr, parameters = c( 'int.lam') )
traceplot( mr, parameters = c( 'beta') )
traceplot( mr, parameters = c( 'eps.i') ) # what parameter is this?

############## whisker plots #############
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#for detection
#fixed effects for detection
whiskerplot( mr, parameters = c( 'int.det', 'alpha' ) , zeroline = TRUE )
#random intercept for observer effect
whiskerplot( mr, parameters = c( "eps.det" ) )
#for abundance
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#fixed effects for abundance
whiskerplot( mr, parameters = c( 'int.lam', 'beta' ) , zeroline = TRUE )
#random intercept for site
whiskerplot( mr, parameters = c( "eps.i" ) )

#derived parameters
whiskerplot( mr, parameters = c( "p" ) )
whiskerplot( mr, parameters = c( "N" ) )



#-----------model evaluation ##########################################

#combine fit statistics into dataframe for plotting
fitdf <- data.frame( fit = mr$sims.list$fit, 
                     fit.new = mr$sims.list$fit.new )
head( fitdf)
#plot
ggplot( fitdf, aes( x = fit, y= fit.new ) ) +
  geom_point( size = 2 ) +
  theme_classic( base_size = 17 ) +
  xlab( "Discrepancy actual data" ) +
  ylab( "Discrepancy predicted data" ) +
  geom_abline(intercept = 0, slope = 1)


#calculate chat:
mean( mr$mean$fit ) / mean( mr$mean$fit.new )



######### violin plots ###########################################
# We may want to plot our model coefficient distributions #
# one approach is to use violin plots
# here we demonstrate for the abundance submodel
#start with extracting relevant fixed effects for abundance submodel
beta.mat <- mr$sims.list$beta
#label columns based on order you included predictors in the model
colnames( beta.mat ) <- c( "elevation", "light" )

# convert matrix to long format
beta.df <- tidyr::gather( data = as.data.frame( beta.mat), key = Predictor,
                          Value )
#plot 
ggplot( beta.df, aes( y = Value, x = Predictor, fill = Predictor ) ) + 
  #set plotting theme 
  theme_classic( base_size = 17 ) + 
  #add label for y axis
  ylab( 'Standardized effect size' ) +
  #plot distribution
  geom_violin( trim = FALSE, fill = '#A4A4A4', 
               size = 0.5 ) +
  #add boxplot highlighting mean and 95 CIs 
  geom_boxplot( width = 0.1 ) +
  #add zero line
  geom_hline( yintercept = 0, size = 1.2, color = 'black' )


alpha.mat<-mr$sims.list$alpha
#colnames( alpha.mat ) <- c( "moon phase" ) #How to name this? there is just one vector

alpha.df<-tidyr::gather(data=as.data.frame(alpha.mat),key=Predictor,Value)

ggplot( alpha.df, aes( y = Value, x = Predictor, fill = Predictor ) ) + 
  #set plotting theme 
  theme_classic( base_size = 17 ) + 
  #add label for y axis
  ylab( 'Standardized effect size' ) + # add xlab 
  #plot distribution
  geom_violin( trim = FALSE, fill = '#A4A4A4', 
               size = 0.5 ) +
  #add boxplot highlighting mean and 95 CIs 
  geom_boxplot( width = 0.1 ) +
  #add zero line
  geom_hline( yintercept = 0, size = 1.2, color = 'black' ) # does this mean the moon phase as a positive effect?


######### partial prediction plots #############################
# Estimate partial prediction plots (marginal effect plots) for predictors 
# with 95% CIs not overlapping zero:
# For abundance submodel first:####
# Start by creating our datasets to predict over
# how many values do we use:
n <- 100# just for the line to be smoth
#define a vector of ones for intercept
int <- rep( 1, n )
# Use the observed values to define range of predictor: # what is this for me?
elevation <- seq( min( s.l.c[,"elevation"]),max(  s.l.c[,"elevation"]),
                  length.out = n )
#standardize predictors:
elevation.std <- scale( elevation ) # we should use the  standardize we use when modeling 

#extract relevant fixed coefficient from abundance sub-model results
fixedabund <- cbind( mr$sims.list$int.lam, mr$sims.list$beta[,2] ) 

#estimate predicted abundance 
predabund <- exp( fixedabund %*% t( cbind( int, elevation.std ) ))
#calculate mean abundance
mabund <- apply( predabund, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
CIabund <- apply( predabund, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
abunddf <- data.frame( mabund, t(CIabund),
                       elevation.std, elevation )
#view
head( abunddf); dim( abunddf)
#rename predicted abundance columns
colnames(abunddf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalised effects for abundance submodel 
ggplot( abunddf, aes( x = elevation, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "Relative abundance" ) +
  xlab( "elevation" ) +
  geom_line( size = 1.5) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI, ymax = highCI ) )



# how to make it for the light variable. 


##### detection marginal effects ######
# our only fixed predictor in detection submodel was time
# what are the min max times:
lano_js %>% select( X25 ,X26, X27,X28, X29, X30, X31, X32, X33 ) %>% #a time for each week?
  summarise_all(list(min, max))

min(obs.cov)
max(obs.cov)

# just do it the same way you did for elevation.


#use them to define your bounds:
Time <- round(seq( 0, 360, length.out = n ),0) ## what should be my bounds?
time.std <- scale( Time )
time2.std <- scale( Time^2 )
#extract relevant fixed coefficient for detection submodel results
fixeddet <- cbind( mr$sims.list$int.det, mr$sims.list$alpha )
#estimate predicted detection
preddet <- plogis( fixeddet %*% t( cbind( int, time.std, time2.std ) ) ) #what is int?
#calculate mean abundance
mdet <- apply( preddet, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
CIdet <- apply( preddet, MARGIN = 2, FUN = quantile, 
                probs = c(0.025, 0.975) )




#create dataframe combining all predicted values for plotting
detdf <- data.frame( mdet, t(CIdet),
                     time.std, time2.std, Time )
#view
head( detdf); dim( detdf)
#rename predicted abundance columns
colnames(detdf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginal effects for abundance submodel 
ggplot( detdf, aes( x = Time, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "Probability of detection" ) +
  xlab( "Time (mins pass 6:00am)" ) +
  geom_line( size = 1.5) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI, ymax = highCI ) )