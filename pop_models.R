

##     created by Carlos Linares            ##
##     
##     Purpose: Replicate analysis from the 
##     https://github.com/quantitativeconservationlab/AppPopnEco/blob/master/CountBayesAnalysis.R
##     
##     and https://github.com/quantitativeconservationlab/AppPopnEco/blob/master/CountBayesPlotting.R
##     

#libraries 
library(tidyverse)
library(jagsUI)

# load data
# this is the data from bat_pop_analysis script and was cleaned in the cleanup_script

js21<-read.csv('data_analysis/bat_pop_analysis/bat_js.csv')
lano_js <-
  read.csv('data_analysis/bat_pop_analysis/lano_js.csv', header = T, col.names = nombres)

#remove the first col from lano because it is not useful. 
lano_js<-lano_js[,-1]
names(lano_js)
nombres <- # was trying to change the name of the columns
  c("",
    "site",
    "25",
    "26",
    "27",
    "28",
    "29",
    "30",
    "31",
    "32",
    "33",
    "34")


#rename week columns as numbers just numbers 



#sites
I<-length(unique(lano_js$site))

#visits
J<-10 # the number of weeks 

# what variables we scale?
# site level covariates

s.l.c<-read.csv('data_analysis/bat_pop_analysis/slc.csv') 
s.l.c<-s.l.c[-1]
# obs covariates

obs.cov<- read.csv('data_analysis/bat_pop_analysis/obs.cov.csv')
obs.cov<-obs.cov[-1]




############## Specify model in bugs language:  #####################
sink( "m1.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
     #random intercept for site
      for( i in 1:I ){
        eps.det[i] ~ dnorm( 0, pres.det ) T(-7, 7)
      }
      #associated variance of random intercepts:     
      pres.det <- 1/ ( sigma.det * sigma.det )
      #sigma prior specified as a student t half-normal:
      sigma.det ~ dt( 0, 2.5, 7 ) T( 0, )
     
       #priors for detection coefficients:
      #define as a slightly informative prior
      #for( a in 1:A){
        alpha ~ dnorm( 0, 0.2 ) T(-7, 7 ) # when we add more obs cov we need to modify and remember to keep one dataframe per covariate
      #}
      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[ b ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #prior for abundance model intercept 
      int.lam ~ dgamma( 0.01, 0.01 ) # this is a vague 
    # ecological model of abundance
    for( i in 1:I ){
      #relative abundance modelled as a Poisson distribution 
      N[ i ] ~ dpois( lambda[ i ] )
      
      #linking abundance rate to predictor
      log( lambda[ i ] ) <- #intercept for abundance 
                          int.lam + 
                 # fixed effects of sagebrush and cheatgrass Fixed for me are vegetation and the light treatment?????
    
                        beta[1]*s.l.c[i,2] + ## this pred is treatment and we add ##elevation and later vegetation add each of the site level pred
                       
    }
      #observation model
     for( i in 1:I ){  
      for( j in 1:J ){
        #model for probability of detection
        logit( p[i,j] ) <- #intercept for detection 
                          int.det + 
                  #random intercept for observer effect
                  eps.det[ i ] + # this makes them random effects.  
    
                  #quadratic effect of time of day # I don't have a quadratic effect of time should we remove this too?????? I just modified it with my own 
    
                  alpha * obs.cov[i,j] 
    
        #observed counts distributed as a Binomial:
        y_obs[ i, j ] ~ dbin( p[i,j], N[i] )  # remember to change it when feeding the data so batjs= y_obs or visceversa. 
                  
        #Model evaluation: We calculate Chi-squared discrepancy
        
        #start with expected abundance
        eval[i,j] <- p[i,j] * N[i]
        #compare vs observed counts
        E[i,j] <- pow( ( y_obs[i,j] - eval[i,j] ), 2 ) /
                  ( eval[i,j] + 0.001 )
        # Generate replicate data and compute fit stats
        #expected counts
        y_hat[ i, j ] ~ dbin( p[i,j], N[i] )
        #compare vs expected counts
        E.new[i,j] <- pow( ( y_hat[i,j] - eval[i,j] ), 2 ) /
                  ( eval[i,j] + 0.001 )
      
    } #close J
    } #close I
        
    #derived estimates of model fit
    fit <- sum( E[,] )
    fit.new <- sum( E.new[,] )
    
    } #model close
     
     ", fill = TRUE )

sink()

################ end of model specification  #####################################     
modelname <- "m1.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.lam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'eps.det' #random intercepts in detection
              , 'sigma.det' #error for random intercept
              , 'beta' #abundance coefficients
              , 'p' #estimate of detection probability
              , 'y_hat' #predicted observations
              , 'N' #estimates of abundance
              , 'fit' #estimate of fit for observed data
              , 'fit.new' #estimate of fit for predicted data
)

#initial values defined as max counts
Nst <- apply( lano_js, 1, max, na.rm = TRUE ) # modify to be batj for one sp

#replace 0s with 1
Nst[which(Nst== 0 )] <- 1

#how many ecological predictors that are fixed effects (s.l.c)
B <- 2
#how many detection predictors that are fixed effects (obs.cov)
A <- 1
#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           N = Nst ) }

#define data that will go in the model
str( win.data <- list( y_obs = as.matrix( lano_js ), # modify to bat x one sp
                       #number of sites, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, #S = S, s=species later 
                       #site level habitat predictors
                       s.l.c = s.l.c,
                       #observation predictors:
                       obs.cov = obs.cov

) )

#call JAGS and summarize posteriors:
m1 <- autojags( win.data, inits = inits, params, modelname, #
                n.chains = 5, n.thin = 10, n.burnin = 20000,
                iter.increment = 10000, max.iter = 500000, 
                Rhat.limit = 1.05,
                save.all.iter = FALSE, parallel = TRUE ) 

#view results 
summary(m1)
plot(m1)
#chat
hist( m1$sims.list$fit / m1$sims.list$fit.new )
#mean chat
mean( m1$mean$fit ) / mean( m1$mean$fit.new )
#Bayesian pvalue
plot( m1$sims.list$fit, m1$sims.list$fit.new )
###### end m1 ########