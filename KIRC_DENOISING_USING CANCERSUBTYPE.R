#main package for this code
library(CancerSubtypes)

library(tidyverse)
library(ggplot2)
library("plyr")

# Preprocess the data
# Load the dplyr package
library("dplyr")

#loading the latent representation dataset & Read the contents of CSV file using read_csv method
KIRC_denoising<- read.csv("latent representations from autoencoder/latentRepresentation_KIRC_Denoising.csv")

#converting into dataframe 

KIRC_denoising<- data.frame(KIRC_denoising)
print(dim(KIRC_denoising))


#replacing all the "TCGA.A2.A04U.01A.11R.A115.07" pattern to "TCGA-A2-A04U-01A-11R-A115-07" for pacients col
KIRC_denoising$Pacients <- gsub(".", "-", KIRC_denoising$Pacients, fixed = TRUE)

# Extract the first 3 parts of each TCGA barcode
KIRC_denoising$Pacients <- substr(KIRC_denoising$Pacients, start = 1, stop = 12)


# write.csv(KIRC_vanilla, "D:/KIRC_VANILLA_PREPROCESSED.csv", row.names=FALSE)



#reading the contents of clinical file
KIRC_clinical_cases<- read_tsv("D:/KIRC_CLINICAL CASES/clinical.tsv",show_col_types = FALSE)

# Convert dataframe to list using data.frame()
KIRC_clinical_cases <- data.frame(KIRC_clinical_cases)


#creating subset of dataframes for the relevant columns
KIRC_filtered_clinical_data<- subset(KIRC_clinical_cases, select=c("case_submitter_id","vital_status","days_to_death","days_to_last_follow_up","ajcc_pathologic_stage"))
KIRC_filtered_clinical_data <- as.data.frame(KIRC_filtered_clinical_data %>% distinct())
print(dim(KIRC_filtered_clinical_data))

###merging "Pacients" and "case_submitter_id" columns from 2 datasets
lr <- KIRC_denoising
cd <- KIRC_filtered_clinical_data
cat(dim(lr), dim(cd))

print(length(intersect(lr$Pacients, cd$case_submitter_id)))
rownames(lr) <- lr$Pacients
rownames(cd) <- cd$case_submitter_id
cat(dim(lr), dim(cd))

print(length(intersect(rownames(lr), rownames(cd))))
merge_data <-merge(lr, cd, by = "row.names" )

print(dim(merge_data))
View(merge_data)

# Removed the repeated columns using select()
merge_data<- merge_data %>% select(-c(Row.names,case_submitter_id))


###filtering the conditions of columns vital_status,days_to_death,days_to_last_follow_up
# merge_data_days_to_death<-filter(merge_data, vital_status == 'Alive' & days_to_death =="'--")
# merge_data_days_to_last_follow_up<-filter(merge_data, vital_status == 'Dead' & days_to_last_follow_up =="'--")

#replacing all the "'--" chars by zero
# merge_data_days_to_death[merge_data_days_to_death == "'--"] <- "0"
# merge_data_days_to_last_follow_up[merge_data_days_to_last_follow_up == "'--"] <- "0"
 

#for a combined df with replaced values
# merge_data <-rbind.fill(merge_data_days_to_death,merge_data_days_to_last_follow_up)
# stages <-merge_data$ajcc_pathologic_stage

###filtering the conditions of columns vital_status,days_to_last_follow_up

# merge_data$vital_status <- as.character(merge_data$vital_status) # convert to character to avoid warnings
# 
# merge_data$days_to_last_follow_up[merge_data$vital_status == "Dead" & is.na(merge_data$days_to_last_follow_up)] <- 0
# 
# merge_data$deceased=merge_data$vital_status == "Dead"
# table(merge_data$deceased)

# create an "overall survival or time" variable that is equal to days to death here
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive


merge_data$time <-ifelse(merge_data$vital_status == "Dead",
                         merge_data$days_to_death,
                         merge_data$days_to_last_follow_up)


merge_data$time = as.numeric(merge_data$time)


### Convert vital_status column to status
merge_data$status <- ifelse(merge_data$vital_status == "Alive", 0, 1)

merge_data<- merge_data %>% select(-c(days_to_death,days_to_last_follow_up))


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
result=ExecuteCC(clusterNum=3,seperate_sub_data_features,maxK=10,clusterAlg="km",distance="euclidean",title="KRCC")
result

group=result$group
distanceMatrix=result$distanceMatrix

#survival analysis
p_value=survAnalysis(mainTitle="KM_Eucli_KRCC_F1000_100F_Denoising_",merge_data$time,merge_data$status,group,
                     distanceMatrix,similarity=TRUE)


#without loop and using PAM and spearman distance

result1=ExecuteCC(clusterNum=3,seperate_sub_data_features,maxK=10,clusterAlg="pam",distance="spearman",title="KRCC")
result1

group=result1$group
distanceMatrix=result1$distanceMatrix

#survival analysis
p_value1=survAnalysis(mainTitle="PAM_spearman_KRCC_F1000_100F_Denoisng",merge_data$time,merge_data$status,group,
                      distanceMatrix,similarity=TRUE)





##in loop
# Create a folder to save the results
dir.create("results_folder")

# Loop through clusterNum values for  k-means and euclidean
for (clusterNum in 3:6) {
  # Execute clustering analysis
  result <- ExecuteCC(clusterNum = clusterNum, seperate_sub_data_features, maxK = 10, clusterAlg = "km", distance = "euclidean", title = "KRCC")
  
  # Extract results
  group <- result$group
  distanceMatrix <- result$distanceMatrix
  
  # Perform survival analysis
  p_value <- survAnalysis(mainTitle = paste0("KM_Eucli_KRCC_F1000_100F_Denoising_"), merge_data$time, merge_data$status, group, distanceMatrix, similarity = TRUE)
  
  # Save the results in a file
  BRCA_Vanilla <- paste0("results_", clusterNum, ".txt")
  filepath <- file.path("results_folder", BRCA_Vanilla)
  write.table(p_value, filepath, row.names = FALSE)
  # ggsave(BRCA_Vanilla,"k means survival analysis on datasets",p_value,device = NULL)
}



# Loop through clusterNum values for PAM and Spearman
for (clusterNum in 3:6) {
  # Execute clustering analysis
  result1 <- ExecuteCC(clusterNum = clusterNum, seperate_sub_data_features, maxK = 10, clusterAlg = "pam", distance = "spearman", title = "KRCC")
  
  # Extract results
  group <- result1$group
  distanceMatrix <- result1$distanceMatrix
  
  # Perform survival analysis
  p_value <- survAnalysis(mainTitle = paste0("PAM_SP_KRCC_F1000_100F_Denoising"), merge_data$time, merge_data$status, group, distanceMatrix, similarity = TRUE)
  
  # Save the results in a file
  BRCA_Vanilla <- paste0("results_", clusterNum, ".txt")
  filepath <- file.path("results_folder", BRCA_Vanilla)
  write.table(p_value, filepath, row.names = FALSE)
  # ggsave(BRCA_Vanilla,"k means survival analysis on datasets",p_value,device = NULL)
}

