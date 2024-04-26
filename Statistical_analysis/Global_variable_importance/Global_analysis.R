###############################################################################################
#
#                                 GLOBAL STATISTICAL ANALYSIS
#
###############################################################################################
# developed by Frederik Saltre 26/04/2024

#This R script is performing a statistical analysis on trajectory data, which includes examining the paths taken by subjects across different environments. Here's a breakdown of the script's key actions in simpler terms:
# Step #1 - Loading Libraries: The script starts by loading two packages: maps, which is used for dealing with geographic data.
#           ade4, which provides tools for multivariate data analysis, often used in ecological studies.
# Step #2 - Reading Data: It reads trajectory data from a CSV file called "AveragedHumanTrajectories(Veghigh)_stats.csv". This data likely contains various measurements for different trajectories taken by subjects, perhaps animals or humans, in a high vegetation area.
# Step #3 - Normalizing Data: A custom function (norme.fct) is defined to normalize the data for each trajectory. Normalizing, in this context, means adjusting the data values so they can be compared more fairly. It adjusts each value based on the average and variability of similar data points (i.e., trajectories that start and end at the same points).
#           This makes it easier to compare trajectories by removing scale differences that come from varying starting and ending points.
# Step #4 - Testing Differences: Another function (test.fct) is set up to test if two groups of trajectories are significantly different from each other. It uses a statistical method called a permutation test, which helps determine if the observed difference between two sets of data could have happened by chance.
#           This is useful for checking, for example, if trajectories that end at a certain point always have a characteristic that is different from those ending elsewhere, which might suggest environmental effects or other influencing factors.
# Step #5 - Loop Through Variables: The script then loops through certain columns of the data (specifically, the 8th to 13th columns). These columns might represent different environmental variables or measurements taken along the trajectories.
#           For each of these variables, it applies the normalization function and then the testing function to see how the trajectory data for the group labeled "Y" compares to other groups.
#           It stores these results in a data frame, which is a type of table in R that makes data easy to manipulate and analyze.
# Step #6 - Output Results: The results, including the variables tested and their corresponding p-values (a statistical measure that helps understand the significance of the results), are printed out.
#           These results help understand which variables significantly differ based on the trajectory group, which could be insightful for ecological or behavioral studies.

#In essence, this script helps analyze how different factors or conditions affect the paths taken by subjects in a high vegetation area, using statistical methods to ensure the findings are robust and not just due to random variations.



#######################
#  LIBRARIES & DATASET
#######################
# Load necessary libraries
library(maps) # For map data
library(ade4) # For ecological and environmental statistical methods in multivariate analysis

# Read in trajectory data from a CSV file
trajm.dat <- read.csv("AveragedHumanTrajectories(Veghigh)_stats.csv",sep=",",dec=".",header=T)

#############
#  FUNCTIONS
#############
# The norme.fct function normalizes the trajectory data based on starting-ending point combinations.
# The test.fct function performs a permutation test to compare the means of two sets of data.

### Function that normalizes the residuals of each trajectory
### by the mean and standard deviation of trajectories starting and ending at the same points
norme.fct <- function(trajm.dat,v0)
{
 v <- NULL
 for(i in unique(trajm.dat[,2]))   # Loop over unique starting-ending combinations
  {
  sel<- trajm.dat[,2]==i # Select rows based on the unique combination
  v1 <- v0[sel]
  # Normalize v1 by subtracting the mean and dividing by the standard deviation
  v1<- (v1-mean(v1))/sd(v1)
  v <- c(v,v1) # Combine results
}
return(v)
}


## Function for permutation test comparing the means of v1 and v2
test.fct <- function(v1,v2)
{
   stat.fct <-function(v1,v2)    # Inner function to calculate difference of means
    { return(mean(v1)-mean(v2))}
   
 obs <- stat.fct(v1,v2) # Observed difference
 sim <- NULL
 # Permutation test simulation
 for(i in 1:10000)
  {
   w <- sample(c(v1,v2)) # Shuffle v1 and v2 together
   sim <- c(sim,stat.fct(w[1:length(v1)],w[(length(v1)+1):length(w)])) # Append simulation results
  }
return(sum(sim>obs)/length(sim))    # Calculate p-value
}


######################
#  MAIN LOOP / PROGRAM
######################
#The main loop iterates over certain columns of the dataset, applying the normalization and testing functions, and prints the results.
# Iterate over environmental variables
# Create a dataframe to store results
out<- data.frame(matrix(ncol = 2, nrow = 0))
colnames(out) <- c('Variable', 'Probablilty')

for (i in 8:13) # Loop through columns 8 to 13 of the data
{
  v <- norme.fct(trajm.dat,trajm.dat[,i]*trajm.dat[,4])   # Normalized values using the norme.fct function
  res <- test.fct(v[trajm.dat[,1]!="Y"],v[trajm.dat[,1]=="Y"])   # Apply the permutation test function
  out<-rbind(out,c(names(trajm.dat)[i],abs(res-0.5)))   # Store results in the dataframe
  print(c(names(trajm.dat)[i],"p-valeur de la diff de moyenne",abs(res-0.5)))
}

# Convert Probability to numeric and scale it
out$Probablilty<-as.numeric(out$Probablilty)
out$scaleprob<-(out$Probablilty - min(out$Probablilty)) / (max(out$Probablilty) - min(out$Probablilty)) 
colnames(out) <- c('Variable', 'Probablilty', 'Scaled Probablilty')
