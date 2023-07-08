#main package for this code
library(CancerSubtypes)
library(tidyverse)
library(ggplot2)
library("plyr")
# Preprocess the data
# Load the dplyr package
library("dplyr")

#set the working directory
#loading the latent representation dataset
latent_rep_data <- read.csv("latent representations from autoencoder/latentRepresentation_GBM_Sparse.csv")
print(dim(latent_rep_data))
#View(latent_rep_data)


# Reading the contents of clinical RDS file using readRDS method
GBM_clinical_data <- readRDS("D:/paper implementation - I/mutliview_subtype_data_code/original_data/GBM_clinical_data.rds")

###merging "Pacients" and "case_submitter_id" 
lr <- latent_rep_data
cd <- GBM_clinical_data
cat(dim(lr), dim(cd))

###setting the index column of clinical data a label
cd <- cbind(case_submitter_id = rownames(cd),cd)
rownames(cd) <- 1:nrow(cd)

#merging on basis of ids
print(length(intersect(lr$Pacients, cd$case_submitter_id)))
rownames(lr) <- lr$Pacients
rownames(cd) <- cd$case_submitter_id
cat(dim(lr), dim(cd))

print(length(intersect(rownames(lr), rownames(cd))))
merge_data <-merge(lr, cd, by = "row.names" )
print(dim(merge_data))
# View(merge_data)

# Remove columns like case submitter id and row names using select()
merge_data <- merge_data %>% select(-c("Row.names","case_submitter_id","race","ethnicity","batch_number","radiations_radiation_regimenindication","dcc_upload_date","karnofsky_performance_score","radiation_therapy","years_to_birth"))

#seperating from merge data the features and the patients ids

seperate_sub_data_features<- merge_data[c(1:101)]


#first column where patient ids were given
#rownames() function used to set the first column as row labels

rownames(seperate_sub_data_features) <- seperate_sub_data_features[, 1]

# Remove the label "patients" from the column names
#colnames(sub_data)[1] <- "Patients"

# Set the first column values to null
seperate_sub_data_features$Pacients <- NULL

# Remove the label "patients" from the column names
#colnames(sub_data)[1] <- ""


#transpose the "patients" column to row labels 
seperate_sub_data_features<- t(seperate_sub_data_features)

#function in R is used to list the available datasets in the loaded packages. 
#It provides a list of datasets that are included with the currently loaded packages.
#data()


#without loop and using k-means and euclidean distance
result=ExecuteCC(clusterNum=3,seperate_sub_data_features,maxK=10,clusterAlg="km",distance="euclidean",title="GBM")
result

group=result$group
distanceMatrix=result$distanceMatrix

#survival analysis
p_value=survAnalysis(mainTitle="GBM_3",merge_data$time,merge_data$status,group,
                     distanceMatrix,similarity=TRUE)


#without loop and using PAM and spearman distance

result1=ExecuteCC(clusterNum=3,seperate_sub_data_features,maxK=10,clusterAlg="pam",distance="spearman",title="GBM")
result1

group=result1$group
distanceMatrix=result1$distanceMatrix

#survival analysis
p_value1=survAnalysis(mainTitle="GBM_3",merge_data$time,merge_data$status,group,
                      distanceMatrix,similarity=TRUE)





##in loop
# Create a folder to save the results
# dir.create("results_folder")

# Loop through clusterNum values for  k-means and euclidean
for (clusterNum in 3:6) {
  # Execute clustering analysis
  result <- ExecuteCC(clusterNum = clusterNum, seperate_sub_data_features, maxK = 10, clusterAlg = "km", distance = "euclidean", title = "GBM")
  
  # Extract results
  group <- result$group
  distanceMatrix <- result$distanceMatrix
  
  # Perform survival analysis
  p_value <- survAnalysis(mainTitle = paste0("GBM_KM_Eucli_F1000_100F_SPARSE"), merge_data$time, merge_data$status, group, distanceMatrix, similarity = TRUE)
  # Save the results in a file
  # BRCA_Vanilla <- paste0("results_", clusterNum, ".txt")
  # filepath <- file.path("results_folder", BRCA_Vanilla)
  # write.table(p_value, filepath, row.names = FALSE)
  # ggsave(BRCA_Vanilla,"k means survival analysis on datasets",p_value,device = NULL)
}



# Loop through clusterNum values for PAM and Spearman
for (clusterNum in 3:6) {
  # Execute clustering analysis
  result1 <- ExecuteCC(clusterNum = clusterNum, seperate_sub_data_features, maxK = 10, clusterAlg = "pam", distance = "spearman", title = "GBM")
  
  # Extract results
  group <- result1$group
  distanceMatrix <- result1$distanceMatrix
  
  # Perform survival analysis
  p_value <- survAnalysis(mainTitle = paste0("GBM_PAM_SPRMN_F1000_100F_Sparse"), merge_data$time, merge_data$status, group, distanceMatrix, similarity = TRUE)
  
  # Save the results in a file
  # BRCA_Vanilla <- paste0("results_", clusterNum, ".txt")
  # filepath <- file.path("results_folder", BRCA_Vanilla)
  # write.table(p_value, filepath, row.names = FALSE)
  # ggsave(BRCA_Vanilla,"k means survival analysis on datasets",p_value,device = NULL)
}
