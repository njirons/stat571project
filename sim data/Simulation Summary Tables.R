## This file is summarize the simulation results for all methods (MLR, MLR+Lasso/Ridge/EN/Entropy)
## This file requires: mae#.RData, mse#.RData,
##                     mape#.RData, mspe#.RData,
##                     scene#list.RData
##                    For Scenario # = 1, 2, 3
##                    scene#list.RData is for MLR, MLR+Lasso/Ridge/EN methods
##                    Other datasets are for MLR+Entropy Methods
## This file outputs: scen#list.RData for scenarios 1,2,3 
##                    list of data.frames with the 4 performance metrics

library(glmnet)
library("MGLM")
library(purrr)
library(nnet)


dir_link <- "/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data/Sim Datasets/"


file_string <- c(sprintf("scen#list.RData",1:3),
                 sprintf("mae%s.RData",1:3),
                 sprintf("mape%s.RData",1:3),
                 sprintf("mse%s.RData",1:3),
                 sprintf("mspe%s.RData",1:3))

## Load simulation datasets
for(i in c(file_string)){
  load(paste0(dir_link,i))
}



###### Compiling MLR + Entropy Statistics #####

# Function to convert all 4 metrics for MLR Entropy into data.frame to match 
# scen#list format for the other methods
ent_metric<-function(x.train, y.train, mse,mae, mspe, mape, nscene){
  #ncat <- ncol(y.train)
  data.frame(dist = nscene,
             sample = nrow(x.train),
             response = ncol(y.train),
             method="Entropy",
             MSE=reduce(mse,c),
             MAE = reduce(mae,c),
             MSPE = reduce(mspe,c),
             MAPE = reduce(mape,c))
}

d<-1
j<-1
ent_metric(xtrain1[[d]][[j]],ytrain1[[d]][[j]],
           mse1[[d]], mae1[[d]],mspe1[[d]],mape1[[d]],1)

## Scenario 1 Entropy
scen1ent_list<- lapply(1:9,   function(d) ent_metric(xtrain1[[d]][[1]],ytrain1[[d]][[1]],
                                                     mse1[[d]], mae1[[d]],mspe1[[d]],mape1[[d]],1)) 
## Scenario 2 Entropy
scen2ent_list<- lapply(1:3,   function(d) cbind(data.frame(rho=rho_vec[d]),ent_metric(xtrain2[[d]][[1]],ytrain2[[d]][[1]],
                                                                                      mse2[[d]], mae2[[d]],mspe2[[d]],mape2[[d]],2)) )
# Scenario 3 Entropy
scen3ent_list<- lapply(1:9,   function(d) ent_metric(xtrain3[[d]][[1]],ytrain3[[d]][[1]],
                                                     mse3[[d]], mae3[[d]],mspe3[[d]],mape3[[d]],3)) 



# Combining Entropy and MLR, MLR+lasso/ridge/EN results together
allscen1_list <- lapply(1:9, function(i) rbind(scen1_list[[i]],scen1ent_list[[i]]))
allscen2_list <- lapply(1:3, function(i) rbind(scen2_list[[i]],scen2ent_list[[i]]))
allscen3_list <- lapply(1:9, function(i) rbind(scen3_list[[i]],scen3ent_list[[i]]))

# Converting data into a data.frame (now each iteration is one row)
allscen1.df <- reduce(allscen1_list, rbind)
allscen2.df <- reduce(allscen2_list, rbind)
allscen3.df <- reduce(allscen3_list, rbind)

########## Summaries ##########

## Using dplyr to get summary statistics

## Scenario 1
scen1summ.df0<- allscen1.df %>% mutate(
  method = factor(method, 
                  levels = c("MLRnnet","Ridge","EN","Lasso","Entropy"),
                  labels =c("MLR","Ridge","EN","Lasso","Entropy")))  %>%
  group_by( sample,response,method) %>% 
  summarise( 
    MSEmean = mean(MSE, na.rm = TRUE),
    MSEsd = sd(MSE, na.rm = TRUE),
    MAEmean = mean(MAE, na.rm = TRUE),
    MAEsd = sd(MAE, na.rm = TRUE),
    MSPEmean = mean(MSPE, na.rm = TRUE),
    MSPEsd = sd(MSPE, na.rm = TRUE),
    MAPEmean = mean(MAPE, na.rm = TRUE),
    MAPEsd = sd(MAPE, na.rm = TRUE)) %>% mutate(
      MSE =paste0(ifelse(MSEmean<0.001, "< 0.001", round(MSEmean,3)), 
                  sprintf(" (%s)", ifelse(MSEsd < 0.001, "< 0.001",round(MSEsd,3) ))),
      MAE =paste0(ifelse(MAEmean<0.001, "< 0.001", round(MAEmean,3)), 
                  sprintf(" (%s)", ifelse(MAEsd < 0.001, "< 0.001",round(MAEsd,3) ))),
      MSPEs =paste0(scientific(MSPEmean, digits=3), 
                    sprintf(" (%s)", scientific(MSPEsd,digits=3) )),
      MSPE =paste0(ifelse(MSPEmean<0.001, "< 0.001", round(MSPEmean,0)), 
                   sprintf(" (%s)", ifelse(MSPEsd < 0.001, "< 0.001",round(MSPEsd,0) ))),
      
      MAPE =paste0(ifelse(MAPEmean<0.001, "< 0.001", round(MAPEmean,0)), 
                   sprintf(" (%s)", ifelse(MAPEsd < 0.001, "< 0.001",round(MAPEsd,0) )))) %>%
  filter(method != "MLRglmnet") #%>% 
#select(-c(MSEmean, MSEsd, MAEmean, MAEsd, MSPEmean, MSPEsd, MAPEmean, MAPEsd))


scen1summ.df <- scen1summ.df0[c("method","sample","response","MAE","MSE","MAPE","MSPE","MSPEs")]

#write.csv(scen1summ.df0, "scenario1full.csv")
#write.csv(scen1summ.df, "scenario 1 summary.csv")

## Scenario 2
scen2summ.df0<- allscen2.df %>% mutate(
  method = factor(method, 
                  levels = c("MLRnnet","Ridge","EN","Lasso","Entropy"),
                  labels =c("MLR","Ridge","EN","Lasso","Entropy")))  %>%
  group_by( sample,response,rho,method) %>% 
  summarise( 
    MSEmean = mean(MSE, na.rm = TRUE),
    MSEsd = sd(MSE, na.rm = TRUE),
    MAEmean = mean(MAE, na.rm = TRUE),
    MAEsd = sd(MAE, na.rm = TRUE),
    MSPEmean = mean(MSPE, na.rm = TRUE),
    MSPEsd = sd(MSPE, na.rm = TRUE),
    MAPEmean = mean(MAPE, na.rm = TRUE),
    MAPEsd = sd(MAPE, na.rm = TRUE)) %>% mutate(
      # method = factor(method, 
      #        levels = c("MLRglmnet","MLRnnet","Ridge","EN","Lasso","Entropy")),
      MSE =paste0(ifelse(MSEmean<0.001, "< 0.001", round(MSEmean,3)), 
                  sprintf(" (%s)", ifelse(MSEsd<0.001, "< 0.001",round(MSEsd,3) ))),
      MAE =paste0(ifelse(MAEmean<0.001, "< 0.001", round(MAEmean,3)), 
                  sprintf(" (%s)", ifelse(MAEsd<0.001, "< 0.001",round(MAEsd,3) ))),
      MSPEs =paste0(scientific(MSPEmean, digits=3), 
                    sprintf(" (%s)", scientific(MSPEsd,digits=3) )),
      MSPE =paste0(ifelse(MSPEmean<0.001, "< 0.001", round(MSPEmean,0)), 
                   sprintf(" (%s)", ifelse(MSPEsd<0.001, "< 0.001",round(MSPEsd,0) ))),
      
      MAPE =paste0(ifelse(MAPEmean<0.001, "< 0.001", round(MAPEmean,0)), 
                   sprintf(" (%s)", ifelse(MAPEsd<0.001, "< 0.001",round(MAPEsd,0) )))) %>%
  filter(method != "MLRglmnet") #%>% 
  #select(-c(MSEmean, MSEsd, MAEmean, MAEsd, MSPEmean, MSPEsd, MAPEmean, MAPEsd))


## including the independent (correlation=0) case to simulation 2 
## with same sample size and response categories
rho0scen1.df <-cbind(data.frame(rho=rep(0,5)),
                     as.data.frame(subset(scen1summ.df0, sample==80&response==50)))[c(2:3,1,4:(ncol(scen1summ.df0)+1))]

scen2summ.df <- rbind(rho0scen1.df,as.data.frame(scen2summ.df0))[c("method","rho","MAE","MSE","MAPE","MSPE","MSPEs")]

#write.csv(scen2summ.df0, "scenario1full.csv")
#write.csv(scen2summ.df, "scenario 1 summary.csv")

## Scenario 3
scen3summ.df0<- allscen3.df %>% mutate(
  method = factor(method, 
                  levels = c("MLRnnet","Ridge","EN","Lasso","Entropy"),
                  labels =c("MLR","Ridge","EN","Lasso","Entropy")))  %>%
  group_by( sample,response,method) %>% 
  summarise( 
    MSEmean = mean(MSE, na.rm = TRUE),
    MSEsd = sd(MSE, na.rm = TRUE),
    MAEmean = mean(MAE, na.rm = TRUE),
    MAEsd = sd(MAE, na.rm = TRUE),
    MSPEmean = mean(MSPE, na.rm = TRUE),
    MSPEsd = sd(MSPE, na.rm = TRUE),
    MAPEmean = mean(MAPE, na.rm = TRUE),
    MAPEsd = sd(MAPE, na.rm = TRUE)) %>% mutate(
      MSE =paste0(ifelse(MSEmean<0.001, "< 0.001", round(MSEmean,3)), 
                  sprintf(" (%s)", ifelse(MSEsd<0.001, "< 0.001",round(MSEsd,3) ))),
      MAE =paste0(ifelse(MAEmean<0.001, "< 0.001", round(MAEmean,3)), 
                  sprintf(" (%s)", ifelse(MAEsd<0.001, "< 0.001",round(MAEsd,3) ))),
      MSPEs =paste0(scientific(MSPEmean, digits=3), 
                    sprintf(" (%s)", scientific(MSPEsd,digits=3) )),
      MSPE =paste0(ifelse(MSPEmean<0.001, "< 0.001", round(MSPEmean,0)), 
                   sprintf(" (%s)", ifelse(MSPEsd<0.001, "< 0.001",round(MSPEsd,0) ))),
      
      MAPE =paste0(ifelse(MAPEmean<0.001, "< 0.001", round(MAPEmean,0)), 
                   sprintf(" (%s)", ifelse(MAPEsd<0.001, "< 0.001",round(MAPEsd,0) )))) %>%
  filter(method != "MLRglmnet") %>% 
  select(-c(MSEmean, MSEsd, MAEmean, MAEsd, MSPEmean, MSPEsd, MAPEmean, MAPEsd))

scen3summ.df <- scen3summ.df0[c("method","sample","response","MAE","MSE","MAPE","MSPE","MSPEs")]

#write.csv(scen3summ.df0, "scenario1full.csv")
#write.csv(scen3summ.df, "scenario 1 summary.csv")

