# stat571project

Data files are in .csv format. 

gvhd_sim.R contains implementations of gradient ascent algorithms for MLE with regular multinomial logistic regression (MLR), as well as entropy-regularized MLR and ridge regression. It also has scripts to load in the GVHD data and fit the models.

fake_data_sim.R generates simulated data and runs the same models on them.
NOTE: the implementation of gradient ascent in fake_data_sim.R is the most up-to-date one and should be used instead of gvhd_sim.R for now
