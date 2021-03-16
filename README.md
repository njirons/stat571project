# stat571project

Data files are in .csv format. 

gvhd_sim.R contains scripts to load in the GVHD data.

fake_data_sim.R generates fake data for the simulation studies. These datasets are stored in the data folder. For example, the data for simulation study 1 are: x1.RData (covariates), y1.RData (response), beta1.RData (coefficients).

ent_mlr_cv.R has functions to fit entropic MLR via fast gradient ascent and to conduct k-fold cross validation.

NOTE: the implementation of gradient ascent in ent_mlr_cv.R is the most up-to-date one and should be used.
