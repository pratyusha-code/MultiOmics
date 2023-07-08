###required libraries
library(CancerSubtypes)
library(dplyr)
library(tidyr)

###read the 3 omicsdatasets

GBM_DNA_data <- readRDS("D:/paper implementation - I/mutliview_subtype_data_code/original_data/GBM_DNA_data.rds")
GBM_miRNA_data <- readRDS("D:/paper implementation - I/mutliview_subtype_data_code/original_data/GBM_miRNA_data.rds")
GBM_mRNA_data <- readRDS("D:/paper implementation - I/mutliview_subtype_data_code/original_data/GBM_mRNA_data.rds")

###convert matrix to data frame
GBM_DNA_data_df <- as.data.frame(GBM_DNA_data)
GBM_mRNA_data_df <- as.data.frame(GBM_mRNA_data)
GBM_miRNA_data_df <- as.data.frame(GBM_miRNA_data)


###transposing columns into rows and rows into columns
GBM_DNA_data_df <-t(GBM_DNA_data_df)
GBM_miRNA_data_df <-t(GBM_miRNA_data_df)
GBM_mRNA_data_df <- t(GBM_mRNA_data_df)

###giving name to the first column of all 3 datasets
GBM_miRNA_data_df <- cbind(Pacients = rownames(GBM_miRNA_data_df),GBM_miRNA_data_df)
rownames(GBM_miRNA_data_df) <- 1:nrow(GBM_miRNA_data_df)

GBM_mRNA_data_df <- cbind(Pacients = rownames(GBM_mRNA_data_df),GBM_mRNA_data_df)
rownames(GBM_mRNA_data_df) <- 1:nrow(GBM_mRNA_data_df)

GBM_DNA_data_df <- cbind(Pacients = rownames(GBM_DNA_data_df),GBM_DNA_data_df)
rownames(GBM_DNA_data_df) <- 1:nrow(GBM_DNA_data_df)

### Extract the first 3 parts of each TCGA barcode of 3 datasets

GBM_miRNA_data_df <- as.data.frame(GBM_miRNA_data_df)
GBM_miRNA_data_df$Pacients<- substr(GBM_miRNA_data_df$Pacients, start = 1, stop = 12)


GBM_DNA_data_df <- as.data.frame(GBM_DNA_data_df)
GBM_DNA_data_df$Pacients <- substr(GBM_DNA_data_df$Pacients, start = 1, stop = 12)

GBM_mRNA_data_df <- as.data.frame(GBM_mRNA_data_df)
GBM_mRNA_data_df$Pacients <- substr(GBM_mRNA_data_df$Pacients, start = 1, stop = 12)


###Identify the common patient IDs in all three datasets using the intersect function
###Find the intersection of patient ids

cat(dim(GBM_DNA_data_df), dim(GBM_miRNA_data_df),dim(GBM_mRNA_data_df))
common_patients <- Reduce(intersect, list(GBM_DNA_data_df$Pacients, GBM_miRNA_data_df$Pacients, GBM_mRNA_data_df$Pacients))
print(length(common_patients))

###miRNA dataset
mirna_subset <-GBM_miRNA_data_df [GBM_miRNA_data_df$Pacients %in%common_patients , ]
print(dim(mirna_subset))
mirna_subset<-t(mirna_subset)
print(dim(mirna_subset))

###mRNA dataset
mrna_subset <- GBM_mRNA_data_df[GBM_mRNA_data_df$Pacients %in%common_patients , ]
print(dim(mrna_subset))
mrna_subset<-t(mrna_subset)
print(dim(mrna_subset))

###dnamethylation_subset
dnamethylation_subset <- GBM_DNA_data_df[GBM_DNA_data_df$Pacients %in%common_patients, ]
print(dim(dnamethylation_subset))
dnamethylation_subset<-t(dnamethylation_subset)
print(dim(dnamethylation_subset))
# cat(dim(mirna_subset), dim(mrna_subset),dim(dnamethylation_subset))




###for miRNA dataset 
column_names <- mirna_subset[1, ]

# Set the first element of column_names to "patients"
# column_names[1] <- "patients"

# Remove the second row
mirna_subset <- mirna_subset[-1, ]

# Set the column names of the matrix
colnames(mirna_subset) <- column_names

# Use fsbyvar function to select features
mirna_select <- FSbyVar(mirna_subset,cut.type="topk", value = 100)
print(dim(mirna_select))


###for mRNA dataset
# colnames(mrna_subset) <- mrna_subset[1, ]

# Extract the first row as column names
column_names <- mrna_subset[1, ]

# Set the first element of column_names to "patients"
# column_names[1] <- "patients"

# Remove the second row
mrna_subset <- mrna_subset[-1, ]

# Set the column names of the matrix
colnames(mrna_subset) <- column_names

# Use fsbyvar function to select features
mrna_select <- FSbyVar(mrna_subset, cut.type="topk",value  = 400)
print(dim(mrna_select))


###for methylation data
# colnames(mrna_subset) <- mrna_subset[1, ]

# Extract the first row as column names
column_names <- dnamethylation_subset[1, ]

# Set the first element of column_names to "patients"
# column_names[1] <- "patients"

# Remove the second row
dnamethylation_subset <- dnamethylation_subset[-1, ]

# Set the column names of the matrix
colnames(dnamethylation_subset) <- column_names

###selecting features using FSbyVar
dnamethylation_select <- FSbyVar(dnamethylation_subset, cut.type="topk", value= 500)
print(dim(dnamethylation_select))

###checking the dimensions
cat(dim(mirna_select), dim(mrna_select),dim(dnamethylation_select))

### Combine the selected features into a single dataframe
mirna_select <- as.data.frame(mirna_select)
mrna_select <- as.data.frame(mrna_select)
dnamethylation_select <- as.data.frame(dnamethylation_select)

merged_data <- bind_rows(mirna_select,mrna_select,dnamethylation_select)
dim(merged_data)

### Transpose the dataframe to get patients as rows and genes as columns
merged_data<-t(merged_data)

### Write the dataframe to a csv file
write.csv(merged_data, "merged_data_GBM.csv")
write.csv(merged_data, "merged_data_GBM1.csv")
