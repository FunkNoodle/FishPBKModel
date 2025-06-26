#           o   o 
#                         /^^^^^7
#           '  '     ,oO))))))))Oo,
#                  ,'))))))))))))))), /{
#             '  ,'o  ))))))))))))))))={
#                >    ))))))))))))))))={        Starting!
#                `,   ))))))\ \)))))))={   
#                  ',))))))))\/)))))' \{
#                    '*O))))))))O*'

### Master Script ###

#This runs each .R script in turn

starttime <- Sys.time()
print(starttime)

# Load the necessary library
library(dplyr)

# Read the PBKinput dataframe
PBKinput <- read.csv("PBKinput.csv")

# Read the template script
template_script <- readLines("PBK Model v6.R")  # Adjust the path to your template

# Iterate over each name in column A
for (name in PBKinput$Name) {
  # Replace the target_name in the template with the current name
  new_script <- gsub('target_name <- "Example"', paste('target_name <- "', name, '"', sep = ""), template_script)
  
  # Create a new filename for the modified script
  new_filename <- paste0(name, ".R")  # Save with the name as part of the filename
  
  # Write the modified script to a new file
  writeLines(new_script, new_filename)
}

source("9Diclofenac.R")
source("10Diclofenac.R")
source("11Diclofenac.R")
source("12Diclofenac.R")
source("13Diclofenac.R")
source("14Diclofenac.R")
source("15Diclofenac.R")
source("16Diclofenac.R")
source("17Diclofenac.R")
source("18Diclofenac.R")
source("19Diclofenac.R")
source("20Diclofenac.R")
source("21Diclofenac.R")
source("22Diclofenac.R")
source("30Fluoxetine.R")
source("31Fluoxetine.R")
source("34Ibuprofen.R")
source("35Ibuprofen.R")
source("36Ibuprofen.R")
source("37Ibuprofen.R")
source("38Ibuprofen.R")
source("39Ibuprofen.R")
source("40Ibuprofen.R")
source("42Ibuprofen.R")
source("43Ibuprofen.R")
source("44Ibuprofen.R")
source("45Ibuprofen.R")

endtime <- Sys.time()
print(endtime)
print(endtime-starttime)

### Writing all blood concentrations to an excel file ###

# List all CSV files in the directory
files <- list.files(pattern = "\\.csv$")

# Create an empty data frame to store results
result <- data.frame(FileName = character(), C_ven = numeric(), C_art = numeric(), stringsAsFactors = FALSE)

# Loop through each file
for (file in files) {
  # Read the CSV file
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Check if the data has rows and the required columns before extracting values
  if (nrow(data) > 0 && all(c("C_ven", "C_art") %in% colnames(data))) {
    # Extract values
    B_value <- data[nrow(data), "C_ven"]
    C_value <- data[nrow(data), "C_art"]
  } else {
    # Set default values if conditions are not met
    B_value <- NA
    C_value <- NA
  }
  
  # Add the file name, C_ven, and C_art values to the result data frame
  result <- rbind(result, data.frame(FileName = file, C_ven = B_value, C_art = C_value, stringsAsFactors = FALSE))
}

# Write the result to a new Excel file using writexl
library(writexl)
write_xlsx(result, "combined_plasma_results.xlsx")

endtime <- Sys.time()

#           o   o 
#                         /^^^^^7
#           '  '     ,oO))))))))Oo,
#                  ,'))))))))))))))), /{
#             '  ,'^  ))))))))))))))))={
#                >    ))))))))))))))))={           Done!
#                `,   ))))))\ \)))))))={   
#                  ',))))))))\/)))))' \{
#                    '*O))))))))O*'

print(endtime)
print(endtime-starttime)