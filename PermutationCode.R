#in this code we simulate the data for the temperature, wind and precipitation 
#in Europe 1000 times with a permutation test using 10000 iterations


#Packages
library(tidyverse) # key package
library(beepr) # Usefuly if you want to add a beep sound and do smth else while the computer is working! :D
library(tseries)# time series multi function package
library(forecast)# same as above
library(ggpubr) #plotting themes
library(timetk) # package for lags
library(zoo)# handles dates and other annoying kind of data
library(lubridate)# extract months and years from dates
library(tsibble)# tidyverse friendly dataframe for time-series
library(mgcv) #GAM models (from leo's master thesis)


setwd("/Users/Soolin/Documents/Niet op desktop/voor nieuwe computer /School/Bachellors thesis/Leonardos code")

#load_data
Dat<-readRDS("dat.Rda")

Out<-readRDS("Out_var.Rda")


##CLEANING AND STANDARIDSATION##

X<-Dat%>%
  mutate(date=as_date(date))%>%
  as_tsibble(index=date)%>%
  fill_gaps()%>%
  tk_augment_lags(c(jet_strength,AO,NAO,US_temp,CAN_temp), .lags = c(1:15), .names = "auto")%>%
  tk_augment_leads(c(jet_strength,AO,NAO,US_temp,CAN_temp), .lags = c(-1:-15), .names = "auto")%>%
  mutate(time_cont=as.integer(as.Date(date, '%Y-%m-%d')-as.Date("1950-01-01")))%>%
  filter(month(date)<3 | month(date)>10)%>%
  mutate(across(starts_with("la"),as.character,.names = "fac_{.col}" ))%>%
  mutate(across(c(where(is.numeric)), scale))%>%
  mutate(across(starts_with("fac"),as.numeric))%>%
  mutate(across(starts_with("fac"),round,digits=2))

rm(Dat)
gc()

Y<-Out%>%
  mutate(date=as_date(date),lon=ifelse(lon==0,360,lon))%>%
  mutate(fac_lat=as.character(lat),fac_lon=as.character(lon))%>%
  as_tsibble(index=date,key=c(lon,lat))%>%
  fill_gaps()%>%
  tk_augment_lags(., c(where(is.numeric)), .lags = c(1:15), .names = "auto")%>%
  tk_augment_leads(., c(4:6), .lags = -1:-15, .names = "auto")%>%
  filter(month(date)<3 | month(date)>10)%>%
  mutate(across(c(where(is.numeric)), scale))%>%
  mutate(fac_lat=as.numeric(fac_lat),fac_lon=as.numeric(fac_lon))%>%
  mutate(across(starts_with("fac"),round,digits=2))
tot<-left_join(Y,X)

gc()

#Simulation function:
# you can choose different effects of the cold spell by changing the values of
# cold_spell_wind_effect,cold_spell_prec_effect,and cold_speel temp_effect
# you can also add different values of autocorrelation for temperatures in North America 
# and all the ouctome variables in Europe.
# The GAM models are there only to maintain some (small) differences between difference locations 
# and a small degree of spatial autocorrelation. 
# Feel free to change different parameters and experiment with it! 


dat_sim<-function(X,Y,seed=23,
                  US_temp_autocorr=0, #fix from the data
                  wind_autocorr=0,
                  prec_autocorr=0.9,
                  temp_autocorr=0.9,
                  cold_spell_wind_effect=1.5,
                  cold_spell_prec_effect=1.5,
                  cold_spell_temp_effect=1.5)
  
  
{
  
  #Gam models to extract small lat,lon effects
  
  wind_mod<-gam(wind~
                  te(lat,lon,
                     k=c(15),
                     d=c(2)),data=tot)
  
  prec_mod<-gam(prec~
                  te(lat,lon,
                     k=c(15),
                     d=c(2)),data=tot)
  
  
  temp_mod<-gam(temp~
                  te(lat,lon,
                     k=c(15),
                     d=c(2)),data=tot)
  
  pred_wind<-predict(wind_mod,type="response")
  pred_prec<-predict(prec_mod,type="response")
  pred_temp<-predict(temp_mod,type="response")
  
  
  #seed for reproductibility
  set.seed(seed)
  
  #Simualtion assuming normality of temperature anomalies, 
  #and weibull distribution with different shape parameters 
  # for wind and precipitation anomalies
  
  
  X_fake<- X%>%
    select(date)%>%
    mutate(US_temp=rnorm(nrow(X)))%>%
    mutate(US_temp=US_temp+US_temp_autocorr*lag(US_temp))%>%
    mutate(US_temp_lag1=lag(US_temp))%>%
    mutate(US_temp_lag2=lag(US_temp_lag1))
  
  Y_fake<- Y%>%
    select(date,lat,lon)%>%
    mutate(wind=pred_wind+rweibull(nrow(Y),shape=3),
           prec=pred_prec+rweibull(nrow(Y),shape=1),
           temp=pred_temp+rnorm(nrow(Y)))
  
  tot_fake<-left_join(Y_fake,X_fake)
  
  #Cold spell effect
  tot_fake<- tot_fake%>%
    mutate(wind=if_else(US_temp_lag2< -2,wind+cold_spell_wind_effect,wind),
           prec=if_else(US_temp_lag2< -2,prec+cold_spell_prec_effect,prec),
           temp=if_else(US_temp_lag2< -2,temp+cold_spell_temp_effect,temp))
  
  tot_fake<-tot_fake%>%
    group_by(lat,lon)%>%
    mutate(wind=wind+wind_autocorr*lag(wind),
           prec=prec+prec_autocorr*lag(prec),
           temp=temp+temp_autocorr*lag(temp))%>%
    ungroup()
  
  return(tot_fake)
  
}

# Input arguments:
# treatment: observations of the treatment, as a vector
# control: observations of contrl, as a vector,
# nsim: number of permutations
mc_permutation_test <- function(treatment, control, nsim){
  
  ## First, get some temp variables that will be used later for speed gain
  n.treat <- length(treatment) # number of observations in the treatment
  n.contr <- length(control) # number of observations in the control
  n <- n.treat + n.contr # total number of observations
  all.data <- c(treatment, control) # treatment and control as one vector
  sum.all <- sum(all.data) # sum of all observations
  
  ## Second, calculate the test statistic for the observed data
  mean_diff <- mean(treatment) - mean(control)
  
  ## Third, we do permutations and compute test statistics.
  mean_perm <- rep(NA, nsim)
  for(i in 1 : nsim){
    
    ## Sample the index to form the permuted treatment group
    treatment_index <- sample(n, size = n.treat, replace = FALSE)
    
    ## treatment and control observations after permutation
    treat_mc <- all.data[treatment_index]
    contr_mc <- all.data[-treatment_index]
    
    ## calculate the difference in mean using the permuted data
    ## we will the identity that sum.all = sum(treat_mc) - sum(contr_mc)
    ## So sum(contr_mc) can be computed much quicker
    sum.treat <- sum(treat_mc) # sum of treatment data
    mean_perm[i] <- sum.treat / n.treat - (sum.all - sum.treat) / n.contr
    
  }
  
  ## Return the observed test statistic and the permuted ones
  list(obs = mean_diff,
       perm = mean_perm)
  
}



Replication <- 1e03 # Number of replications
N.Perm <- 1e04 # Number of permutations
## Some variables from your code
locations <- 1 #What locations to look at
colnr_var_out <- 3 #1=temp 2=wind 3=prec
colnr_var_in <- 4 #US temperature
colnr_locations=5 #vector containing location data
colnr_index=9 #vector containing indexes
alpha=0.05 
lag=5 #lag from a cold spell in the US to the effect in Europe 

#Define a empty vector for later
Proportion <- rep(NA, Replication)
for(reps in 1 : Replication){
  
  #generate the data
  tot_fake<-dat_sim(X,Y,seed=23 * reps)
  
  dat_fake_no_effect<-tot_fake%>%
    select(date,US_temp)%>%
    distinct(date,.keep_all=TRUE)
  
  out_fake_no_effect<-tot_fake%>%
    select(date,lat,lon,wind,prec,temp)
  
  
  ################Formatting data####################
  dat <- dat_fake_no_effect
  #sorting the data by temperature 
  #first add an index corresponding dates
  dat<-dat%>%
    group_by(date)%>%
    mutate(index=cur_group_id())%>%
    ungroup()
  
  dat_sort <- dat[order(dat$US_temp),]
  dat_sort$index <- dat$index[order(dat$US_temp)] #working with indexes is easier 
  #as the as.Date() function brings many complications
  
  i <- 1 #starting position for the index
  x_temp <- dat_sort$US_temp[1] #start for coldest temp
  x_date <- dat_sort$index[1] #starting position for coldest date
  x_date_to_remove <- c() #just giving it a starting value
  
  while(i < 0.05*length(dat_sort$index)){ #choosing the 5% coldest days 
    #checking if a cold spell occurred within +-5 days 
    if( (dat_sort$index[i]-5) %in% x_date | (dat_sort$index[i]-4) %in% x_date | (dat_sort$index[i]-3) %in% x_date
        | (dat_sort$index[i]-2) %in% x_date | (dat_sort$index[i]-1) %in% x_date | (dat_sort$index[i]) %in% x_date
        | (dat_sort$index[i]+1) %in% x_date | (dat_sort$index[i]+2) %in% x_date | (dat_sort$index[i]+3) %in% x_date
        | (dat_sort$index[i]+4) %in% x_date | (dat_sort$index[i]+5) %in% x_date){ #not counting cold spells within +-5 days 
      i <- i+1
      #saving the dates that occur with 5 days of another cold spell, to remove them later
      x_date_to_remove <- c(x_date_to_remove, dat_sort$index[i])
    }
    else{
      #adding cold spells only if they occur later than 5 days after or before the nearest cold spell
      x_temp <- c(x_temp, dat_sort$US_temp[i])
      x_date <- c(x_date, dat_sort$index[i])
      i <- i+1
    }
  }
  
  x_temp #all cold spell temperatures
  x_date #all cold spell dates
  x_date_to_remove #dates that occur within 5 days of a cold spell
  
  dat <-  dat %>% filter(!dat$index %in% x_date_to_remove) #removing the observations where cold days
  #occur within 5 days of another cold day
  
  out <- out_fake_no_effect
  
  #converting character dates to numeric values
  out$date <- as.Date(out$date)
  out$day <- day(out$date)
  out$month <- month(out$date)
  out$year <- year(out$date)
  
  #creating a column with all unique dates going from 1 to n
  out<-out%>%
    group_by(date)%>%
    mutate(index=cur_group_id())%>%
    ungroup()
  
  #removing the observations where a 5% coldest day was within 5 days of another cold day
  out <- out %>% filter(!out$index %in% (x_date_to_remove+lag))
  
  #putting everything into one data frame by adding each US temperature once for every latlon combination
  out<-out%>%
    left_join(dat)
  
  #grouping latitute and longitute by location
  out<-out%>%
    group_by(lat,lon)%>%
    mutate(location=cur_group_id())%>%
    ungroup()
  
  #removing october and march
  out <- out[out$month<3 | out$month>10 ,]
  
  #centering the variables around zero as they are not
  #Standardize the variables based on location as each location has different variance
  #create a subset holding the variable we want to scale and the location
  
  out<-out%>%
    group_by(location)%>%
    mutate(temp=scale(temp),wind=scale(wind),prec=scale(prec))%>%
    ungroup()
  
  #creating a matrix as matrices are faster than data frames
  data <- as.matrix(cbind(out$temp,
                          out$wind,
                          out$prec,
                          out$US_temp,
                          out$location,
                          out$day,
                          out$month,
                          out$year,
                          out$index))
  colnames(data) <- c("temp", "wind", "prec", "US_temp",
                      "location", "day", "month", "year", "index")
  
  #removing some dates where the simulation has NA values
  data <- na.omit(data)
  print(data[1: 2,])
  
  ##------------------------------------------------------------------------##
  ## Perform permutation test
  ##------------------------------------------------------------------------##
  ## First, extract the temperature after cold spell
  Treat.Index <- which(data[,colnr_index] %in% (x_date+lag) 
                       & data[,colnr_locations] %in% locations)
  treatment <- data[Treat.Index, colnr_var_out]
  ## Second, extract the temperature that does not correspond to the one after cold spell
  Contr.Index <- which(!data[,colnr_index] %in% (x_date+lag) 
                       & data[,colnr_locations] %in% locations)
  control <- data[Contr.Index, colnr_var_out]
  ## Third, we perform the permutation test
  Perm <- mc_permutation_test(treatment = treatment, control = control, nsim = N.Perm)
  ## Forth, decide whether to reject or not
  
  #two sided test
  #Proportion[reps] <- mean(abs(Perm$perm) > abs(Perm$obs))
  #one sided test
  Proportion[reps] <- mean(Perm$perm > Perm$obs)
  
  
  ## Check simulation progress
  cat("Replication =", reps, "\n")
  
}

hist(Proportion)
mean(Proportion <= 0.05) #How many times the test was rejected
beep(5)

