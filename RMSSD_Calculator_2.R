# 
# |**********************************************************************;
#   * Project           : SRMSSD_Calculator
# *
#   * Program name      : SRMSSD_Calculator.R
# *
#   * Author            : Alex Pavlides
# *
#   * Date created      : 12/08/2016
# *
#   * Purpose           : Calculate RMSSD from IBT values, for time windows related to stimulus events of participants.
#                         Output includes a table of RMSSD value per participant and per stimulus event called RMSSD_Table.csv
#                         Output also includes a figure showing the mean RMSSD over participants for events greater than 15 seconds. 
# *
#   |**********************************************************************;
# *   Parameters        :
# *   path - path to data folder                                                     
# *   num_cores - number of cores to use in parallelised code                                                    
# * 
# *   Example call      :                                                                          
# *   setwd("D:/APavlides/Documents/Kaggle/GitHub Data") -set working directory that contains function
# *   source("RMSSD_Calculator_2.R") - load function  
# *   RMSSD_Calculator_2("D:/APavlides/Documents/Kaggle/GitHub Data/join_us-correction/data_engineer/data", 7) -- run function

RMSSD_Calculator_2 <- function(path, num_cores){

#Packages
library('plyr')
library('foreach')
library('doSNOW')  

cl<-makeCluster(num_cores) #change the 2 to your number of CPU cores  
registerDoSNOW(cl)  

#Retrieve paths for data and the ID's of participants
files   <- list.files(path,recursive = TRUE,full.names=F)
files   <- files[order(as.numeric(sub("([0-9]*).*", "\\1", files)))]
file_ID <- list.files(path,full.names=F)
file_ID <-file_ID[order(as.numeric(sub("([0-9]*).*", "\\1", file_ID)))]

#init vectors 
file_ID_4Deletion <- vector(mode="integer", length=0)
stim_time_window <- vector(mode="numeric", length=0)
RMSSD <- vector(mode="numeric", length=0)
names_header <- vector(mode="character", length=0)
names_header_filtered <- vector(mode="character", length=0)

#initialise list
RMSSD_Store <- list()
RMSSD_Store_filtered <- list()

count_ID <- 1
  #remove incomplete data where there is not both the sensor AND stimulation data.
  for (i in seq(1,length(file_ID),1)){
    
    FindString <- paste(paste ("^",file_ID[i], sep ="", collapse = NULL), "/", sep="", collapse = NULL)
    occurrence <- length(grep(FindString, files))
  
    if (occurrence < 2){
    ToRemove <- grep(FindString, files)
    files <-files[-ToRemove]
    file_ID_4Deletion[count_ID] <- as.integer(file_ID[i])
    count_ID = count_ID + 1
    }
  }

#remove subject ID's with incomplete data
file_ID <- file_ID[! file_ID %in% file_ID_4Deletion]

#init number of participants
Number_of_participants <- (length(files)/2)

comb <- function(x, ...) {
  lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#main loop for calculating RMSSD from IBT
output <- foreach(m = 1:Number_of_participants,.combine='comb', .multicombine=TRUE,
        .init=list(list(), list())) %dopar% {
          
          #retrieve data for subject
          hr_gsr_path   <- paste (path,files[2*m-1], sep ="/", collapse = NULL)
          stimulus_path <- paste (path,files[2*m], sep ="/", collapse = NULL)
          hr_gsr   = read.csv( hr_gsr_path, header = TRUE)
          stimulus = read.csv( stimulus_path, header = TRUE)
          
          #convert time formats
          #convert unix utc to time, using Central European Time Zone
          hr_gsr$unix_utc <- as.POSIXct(as.numeric(as.character(hr_gsr$unix_utc)),origin="1970-01-01",tz="GMT")
          names(hr_gsr)[names(hr_gsr) == 'unix_utc'] <- 'sensor_time'
          
          #remove atypical IBT
          lowerNormal <- which(hr_gsr[,3] < 600) 
          higherNormal <- which(hr_gsr[,3] >1200)
          removeAtypical <- c(lowerNormal,higherNormal)
          # catch integer(0) error
          if (length(removeAtypical) != 0){ 
          hr_gsr <- hr_gsr[-removeAtypical,]}
                 
          #change stimulus time format
          stimulus$server_time <- as.POSIXct(stimulus$server_time,origin="1970-01-01",tz="GMT")
          names(stimulus)[names(stimulus) == 'server_time'] <- 'stimulus_time'
          
          #find no. of rows
          Nrows_hr_gsr <- nrow(hr_gsr)
          Nrows_stimulus <- nrow(stimulus)
          
          stim_count <- 1
          #Get time window length (sec) of stimulus
          for (i in seq(1,Nrows_stimulus,2)){      
            stim_time_window[stim_count] <- as.numeric(stimulus[i+1,1] - stimulus[i,1])
            stim_count = stim_count + 1
          }
          
          #Calculate RMSSD for each stimulus window 
          count_RMSSD <- 1  
          for (i in seq(1,Nrows_stimulus,2)){
            
            start_time <- stimulus[i,1]
            end_time   <- stimulus[i+1,1]              
            RMSSD_sum <- 0
            
            for (n in seq(1,Nrows_hr_gsr-1,1)){
              
              #if sensor time is between stimulation times calculate MRSSD_sum
              if (hr_gsr[n,1] >= start_time && hr_gsr[n,1] < end_time){RMSSD_sum <- RMSSD_sum  +  (hr_gsr[n+1,3] - hr_gsr[n,3])^2}
              #if sensor time is less than start of stimulation time, continue to next time step
              else if(hr_gsr[n,1] < start_time){next} 
              #if sensor time is greater than stimulation time, break and calculate RMSSD
              else break 
            }
            
            #calculate the RMSSD using RMSSD_sum
            RMSSD[count_RMSSD] <- sqrt(RMSSD_sum*(1/(Nrows_hr_gsr-1)))
            
            #reset sum 
            RMSSD_sum <- 0
            
            #iterate count
            count_RMSSD <- count_RMSSD + 1
          }
          
          #list of lists (foreach can only output a list so bundle two lists together for output)
          output <- list(RMSSD,  RMSSD[stim_time_window>15])
          
} # end of foreach
  
stopCluster(cl)  
  
  #extract first list, RMSSD
  RMSSD_Store <- output[[1]]
  ## Add NA values to list elements
  l1 <- lapply(RMSSD_Store , `length<-`, max(lengths(RMSSD_Store )))
  output_RMSSD_Store <- as.data.frame(l1,stringsAsFactors=FALSE)
  
  #extract second list, RMSSD where stimulation events longer than 15 seconds
  RMSSD_Store_filtered <- output[[2]]
  ## Add NA values to list elements
  l2 <- lapply(RMSSD_Store_filtered , `length<-`, max(lengths(RMSSD_Store_filtered )))
  output_RMSSD_Store_filtered <- as.data.frame(l2,stringsAsFactors=FALSE)
      
      #create headers for data frame
      for (k in seq(1,Number_of_participants,1)){
          names_header[k] <- paste ('RMSSD', file_ID[k], sep =":", collapse = NULL)            
      }
      
      #add column headers to data frame
      names(output_RMSSD_Store) <- names_header
      names(output_RMSSD_Store_filtered) <- names_header
      
      #calculate mean RMSSD by subject for each stimulation event, avoiding NA by setting na.rm = T      
      output_RMSSD_Store_filtered[output_RMSSD_Store_filtered == 0] <- NA #remove what is most likely spurious zeros in the data.
      output_RMSSD_Store_filtered <- output_RMSSD_Store_filtered[rowSums(is.na(output_RMSSD_Store_filtered))!=ncol(output_RMSSD_Store_filtered), ] #trim all NA rows
      stim_mean <- vector(mode="numeric", length=0)
      for (n in 1:nrow(output_RMSSD_Store_filtered)){     
        stim_mean[n] <- rowMeans(output_RMSSD_Store_filtered[n,], na.rm = TRUE)
      }

  #output to .csv
  write.csv(output_RMSSD_Store, file = "RMSSD_Table.csv")
  
  #create line chart
  OutputRowN <- nrow(output_RMSSD_Store_filtered)
  stim_events <- seq(1,OutputRowN,1)
  xrange <- range(stim_events)
  yrange <- range(stim_mean)
  plot(xrange, yrange, type="n", xlab="Stimulation Event", ylab="Mean RMSSD of Event" )
  lines(stim_events, stim_mean, type="b", lwd=1.5)
  title("Mean RMSSD across participants associated with Stimulation Events > 15 sec", cex.lab=0.75)
  
}


 