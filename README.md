# SynetiqTask
This function calculate RMSSD from IBT values, for time windows related to stimulus events of participants. The output includes a table of RMSSD value per participant and per stimulus event called RMSSD_Table.csv. In addition, the output also includes a figure showing the mean RMSSD over participants for events greater than 15 seconds.

The function(path, num_cores) parameters are:
*   path - path to data folder                                                     
*   num_cores - number of cores to use in parallelised code

Example call:                                                                          
*   setwd("D:/APavlides/Documents/Kaggle/GitHub Data") -set working directory that contains function
*   source("RMSSD_Calculator_2.R") - load function  
*   RMSSD_Calculator_2("D:/APavlides/Documents/Kaggle/GitHub Data/join_us-correction/data_engineer/data", 7) -- run function



