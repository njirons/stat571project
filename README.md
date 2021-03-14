# stat571project

Data files are in .csv format. 

gvhd_sim.R contains implementations of gradient ascent algorithms for MLE with regular multinomial logistic regression (MLR), as well as entropy-regularized MLR and ridge regression. It also has scripts to load in the GVHD data and fit the models.

fake_data_sim.R generates simulated data and runs the same models on them.

ent_mlr_cv.R has functions to fit entropic MLR via fast gradient ascent and to conduct k-fold cross validation.

NOTE: the implementation of gradient ascent in ent_mlr_cv.R is the most up-to-date one and should be used.
