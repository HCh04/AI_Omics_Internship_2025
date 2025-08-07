# ===================================================================
#             AI and Omics Resrarch Internship (2025)
#             Module I: Getting Started with R (Class Ib)
#            #### notes made for me to follow easier ####
# ===================================================================

#### Set up folders (point of clarity) ####
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

#### Load Dataset ####
data<-read.csv("raw_data/patient_info.csv")

#### Inspect Dataset ####
str(data) # This shows structure: var names + data types
summary(data) # This shows descriptive stats for each var

#### Fix Inconsistent Data Types ####

# Convert 'gender' to factor
data$gender <- as.factor(data$gender)
data$gender <- ifelse(data$gender == "Female", 1, 0)
data$gender <- as.factor(data$gender)

# Convert diagnosis to factor
data$diagnosis_fac <- as.factor(data$diagnosis)
data$diagnosis_fac <- factor(data$diagnosis_fac, levels = c("Normal", "Cancer"))

# Convert age and bmi to numeric
data$age <- as.numeric(data$age)
data$bmi <- as.numeric(data$bmi)

# Convert 'smoker' to factor
data$smoker_fac <- as.factor(data$smoker)
#### Create Binary Smoking Status ####
data$smoker_binary <- ifelse(data$smoker == "Yes", 1, 0)
data$smoker_binary <- as.factor(data$smoker_binary)

#### Saving Cleaned Data ####
dir.create("clean_data", showWarnings = FALSE)
write.csv(data, "clean_data/patient_info_clean.csv", row.names = FALSE)

#### Save workspace ####
save.image(file = "results/HasmikChilingaryan_ClassIb_Assignment.RData")
