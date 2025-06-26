#           o   o 
#                         /^^^^^7 
#           '  '     ,oO))))))))Oo,
#                  ,'))))))))))))))), /{ ~~~~~~~~~~
#             '  ,'o  ))))))))))))))))={  ~~~~~~~~~~   Speedyfish
#                >    ))))))))))))))))={    ~~~~~~~~~~   Starting!
#                `,   ))))))\ \)))))))={  ~~~~~~~~~~
#                  ',))))))))\/)))))' \{ ~~~~~~~~~~
#                    '*O))))))))O*' 

# This runs each .R script in parallel
# Faster than the Batch Code, but will slow your computer down

### Master Script ###

starttime <- Sys.time()
print(starttime)

# Load necessary libraries
install.packages(c("doParallel", "foreach"))
library(doParallel)
library(foreach)
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

# List of scripts to source
scripts <- c(
  
  "6Carbamazepine.R",
  "7Carbamazepine.R",
  "8Carbamazepine.R",
  "11Diclofenac.R",
  "12Diclofenac.R",
  "13Diclofenac.R",
  "14Diclofenac.R",
  "15Diclofenac.R",
  "19Diclofenac.R",
  "20Diclofenac.R",
  "21Diclofenac.R",
  "22Diclofenac.R",
  "29Diphenhydramine.R",
 "34Ibuprofen.R",
  "35Ibuprofen.R",
  "36Ibuprofen.R",
  "37Ibuprofen.R",
  "38Ibuprofen.R",
  "39Ibuprofen.R",
  "40Ibuprofen.R"
  

)

# Set up parallel backend
num_cores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run the scripts in parallel
foreach(script = scripts) %dopar% {
  source(script)
}

# Stop the cluster
stopCluster(cl)


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