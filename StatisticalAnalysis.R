
# This R Script is the Statistical Analysis carried out by the Statistics
# subgroup. The aim is to look at the Geographic comparison between whistle
# categories in Brazil, the Mediterranean and Hawaii. The Sotalia was used
# to guide the methods except in this case our species is equivalent to whistle category. 
# Any questions please email Sophie Gibson at sjg27@ or Danielle Harris at dh17@.



# Packages needed after reading the Sotalia paper
library(ade4)
library(divo)
library(abdiv)
library(iNEXT)
library(SpadeR)
library(vegan)
library(geosphere)
library(ape)

####################          Data Manipulation          ######################

# Create a data frame with two columns, pooled data name and whistle category 
# Order the categories from 1 to 94.

df <- data.frame(pooled_data$name,pooled_data$category) # create data frame
ordered_df <- df[order(df$pooled_data.category),] # order categories

# Create a data frame which is Type(1) abundance data for the SimilarityMult 
# function which includes the three locations and the number of whistles in 
# categories 1 to 94 for each location. 
# This was calculated manually by using ordered_df and counting the number of 
# whistles in each category for each location (stored in Excel sheet too)

# Number of whistles in order of category 1 to 94 for Brazil

brazil_whistles <- c(5,4,1,0,0,9,1,0,0,0,0,5,12,8,0,0,1,0,0,0,2,0,0,0,0,0,2,0,0
                     ,0,0,0,0,4,0,0,0,2,0,0,0,1,2,0,1,4,5,0,0,0,0,0,0,0,7,5,0,
                     0,0,1,1,0,0,3,0,0,1,4,0,0,0,0,0,0,0,1,0,4,0,2,0,0,1,1,0,0,0
                     ,0,0,0,0,2,1,4)

# Number of whistles in order of category 1 to 94 for Med

med_whistles <- c(0,3,0,0,0,1,7,3,0,12,5,8,3,0,0,1,0,4,1,3,1,0,0,1,4,0,1,0,8,1
                  ,0,0,0,0,8,2,0,0,0,0,1,0,0,1,1,3,0,1,0,0,0,2,0,0,1,1,0,0,0,0,
                  0,0,0,2,0,0,0,3,2,1,2,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,4,0,0)

# Number of whistles in order of category 1 to 94 for Hawaii

hawaii_whistles <- c(2,1,0,1,1,4,3,0,2,0,0,0,2,4,2,4,0,0,0,0,0,1,1,10,1,1,1,1
                     ,0,0,1,1,1,2,0,1,2,0,2,2,0,2,0,1,2,2,1,0,1,2,1,0,4,1,2,1
                     ,1,2,1,0,2,1,1,1,1,2,0,2,0,0,0,1,0,1,1,0,1,3,1,1,1,2,0
                     ,0,1,1,1,1,1,1,1,0,0,0)

# Create a data frame for the SimilarityMult function which includes the number of 
# whistles in each category for Brazil,the Mediterranean and Hawaii.

whistles_df <- data.frame(brazil_whistles,med_whistles,hawaii_whistles)
whistles_df

# Estimate various 3 location similarity indices. The abundance- based indices
# include the Horn and Morisita Horn indices.

help("SimilarityMult")
help("Genetics")

# I found the Genetics output gives the 1- pairwise similarity matrix that is 
# found in SimilarityMult (i.e. pairwise dis-similarity matrix). Although, still
# used SimilarityMult to gain Horn and Morisita Horn estimates.

SimilarityMult(whistles_df, datatype='abundance') # look at SimilarityMult output
test_h <- Genetics(whistles_df,q=1) # get Horn matrix
test_mh <- Genetics(whistles_df,q=2) # get Morisita Horn matrix


# From the SimilarityMult output the Horn estimate is 0.4825 +/- 0.1593 and 
# Morisita-Horn is 0.3673 +/- 0.0875.
# The lower similarity values at higher orders of q suggests increased differentiation
# in more common whistle categories.


# Extract a pairwise dis-similarity matrix from the Genetics output for Horn
# and Morisita Horn.

pairwise_h <- test_h$dissimilarity_matrix$Horn # extract Horn dis-similarity matrix
pairwise_mh <- test_mh$dissimilarity_matrix$C22 # extract Morisita Horn dis-similarity matrix
class(pairwise_h) # check it is a matrix
class(pairwise_mh) # check it is a matrix

############ Visualizing pairwise similarity between populations ##############

# The Sotalia paper used the function metaMDS but it gives our data set a warning message
# and a stress of zero which suggests our sample size is too small hence
# a recommendation is to use metric scaling like cmdscale and wcmdscale.
# The only difference between them is one is weighted and one is not.
# I get the same answer for both as our sample is equally weighted.


# Look at vegan function and help files for the scaling functions.
help(vegan)
help("metaMDS")
help("cmdscale")
help("wcmdscale")

# Using Weighted Classical(Metric) Scaling for the Horn Plot
sol_w_horn<- wcmdscale(pairwise_h)
sol_w_horn

# Compare to Classical(Metric) Multidimensional Scaling. Since it gives the
# same result, choose either methods.
sol_c_horn<- cmdscale(pairwise_h)
sol_c_horn

# Using Weighted Classical(Metric) Scaling For Morisita Horn Plot
sol_w_mhorn<- wcmdscale(pairwise_mh)
sol_w_mhorn


####################       Creating Hull Polygon      #########################

# Overlay hull polygons

Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

###############   Plotting the Horn and Morisita Horn estimates  ##############

# Plot Horn estimates
plot_w_horn <- plot(sol_w_horn)
Plot_ConvexHull(sol_w_horn[,1],sol_w_horn[,2],lcolor="yellow") # For Morisita Horn plot

# Plot Morisita Horn estimates
plot_w_mhorn <- plot(sol_w_mhorn)
Plot_ConvexHull(sol_w_mhorn[,1],sol_w_mhorn[,2],lcolor="magenta") # For Morisita Horn plot

# Explore the ggplot library to replicate the plots in the Sotalia paper.
library(ggplot2)


######################     Mantel tests      #########################

# Investigate the inputs needed for the Mantel test.
help("mantel")

# Geographic distance calculation using the Haversine method
help("distHaversine")

# Use the distHaversine function to get the distance between our three 
# geographical points according to the Haversine formula.
# Create a distance matrix for this, however since the data collection locations
# for each location are spread out, we need to decide what coordinate to use.
# Therefore, the mantel test will be disregarded until further information on
# data collection locations are found. Once this happens call the following
# function: mantel(pairwise_h, future_distance_matrix,method ='spearman') when
# using the Horn dis-similairty matrix.




#####################         Regression Analysis         #####################

# Explore the iNEXT package.
?iNEXT
# Look at examples given in the help file
data(spider)
out1 <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
out1$DataInfo
out1$AsyEst$Estimator

# Now using the whistles data frame from earlier, find estimates to use in 
# the regression analysis for each location.

iNEXT_output <- iNEXT(whistles_df, q=c(0,1,2), datatype="abundance")

# Plot Species diversity vs Number of individuals
ggiNEXT(iNEXT_output, facet.var="Assemblage", color.var="Order.q") # q values on same graph
ggiNEXT(iNEXT_output, facet.var="Both", color.var="Order.q") # seperate graph for each q value

# Extract relevant indices from iNEXT output. This output contains values for
# Species richness, Shannon and Simpson indices (q = 0, 1 and 2 respectively)
iNEXT_output$AsyEst 

# However, the Statistics subgroup have decided we can not perform a sensible
# regression with 3 data points. If there were more locations, that could be an 
# option.

