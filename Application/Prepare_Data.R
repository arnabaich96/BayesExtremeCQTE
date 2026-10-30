# Downloads the Vanderbilt diabetes data, defines obesity (BMI >= 30) and imputes missing covariates;
# BMI defines the treatment only and is not kept. Output: Application/data/Diabetes_Data.rds

# --- 1. Load All Required Packages ---
suppressPackageStartupMessages({
  library(here)
  library(mice)
})


# --- 2. Diabetes Data Preparation ---
print("--- Starting Diabetes Data Preparation ---")

# Raw data, downloaded once from the Vanderbilt Biostatistics Datasets
data_dir <- here("Application/data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
raw_file <- file.path(data_dir, "diabetes.csv")
if (!file.exists(raw_file)) {
  download.file("https://hbiostat.org/data/repo/diabetes.csv", raw_file, mode = "wb")
}
diabetes_raw_data <- read.csv(raw_file)
# ratio and glyhb rounded to 10 significant digits, as in the copy of the data used for the article
sig10 <- function(x) {
  ok <- !is.na(x)
  x[ok] <- as.numeric(sprintf("%.10g", x[ok] * (1 + 1e-14)))
  x
}
diabetes_raw_data$ratio <- sig10(diabetes_raw_data$ratio)
diabetes_raw_data$glyhb <- sig10(diabetes_raw_data$glyhb)
cat("Raw data:", nrow(diabetes_raw_data), "subjects\n")

# Outcome and treatment must be observed
diabetes_raw_data <- subset(diabetes_raw_data, !is.na(glyhb) & !is.na(height) & !is.na(weight))
cat("With observed HbA1c, height and weight:", nrow(diabetes_raw_data), "subjects\n")

# Convert categorical columns to numeric codes
diabetes_raw_data$gender   <- as.numeric(factor(diabetes_raw_data$gender))                 # 1 = female, 2 = male
diabetes_raw_data$frame[diabetes_raw_data$frame == ""] <- NA                               # blank = missing
diabetes_raw_data$frame    <- as.numeric(factor(diabetes_raw_data$frame,
                                                levels = c("small", "medium", "large")))   # ordinal 1-3
diabetes_raw_data$location <- as.numeric(factor(diabetes_raw_data$location))               # 1 = Buckingham, 2 = Louisa

# Treatment
bmi <- 703 * diabetes_raw_data$weight / diabetes_raw_data$height^2
diabetes_raw_data$obese <- as.integer(bmi >= 30)

# MICE with predictive mean matching; the analysis uses the first completed data set (m = 10 as in the article's run)
analysis_vars <- c("glyhb", "obese", "stab.glu", "ratio", "location", "age", "gender",
                   "frame", "bp.1s", "bp.1d", "time.ppn")
print("Running MICE imputation for Diabetes data... (This may take a moment)")
print(colSums(is.na(diabetes_raw_data[, analysis_vars])))
imputed_data <- mice::mice(diabetes_raw_data[, analysis_vars], m = 10, maxit = 20, method = "pmm",
                           seed = 123, printFlag = FALSE)
Diabetes_Data_Clean <- cbind(id = diabetes_raw_data$id, complete(imputed_data, 1))

# --- Print Cleaned Data before Saving ---
print("--- Displaying top 6 rows of the cleaned Diabetes data ---")
print(head(Diabetes_Data_Clean))
cat("\n")
cat("n =", nrow(Diabetes_Data_Clean), "| obese (BMI >= 30):", sum(Diabetes_Data_Clean$obese),
    "| HbA1c > 7:", round(100 * mean(Diabetes_Data_Clean$glyhb > 7), 1), "%\n")
cat("BMI range: non-obese", round(range(bmi[diabetes_raw_data$obese == 0]), 2),
    "| obese", round(range(bmi[diabetes_raw_data$obese == 1]), 2), "\n")
print(table(gender = Diabetes_Data_Clean$gender, obese = Diabetes_Data_Clean$obese))


saveRDS(Diabetes_Data_Clean, file.path(data_dir, "Diabetes_Data.rds"))

print("--- Diabetes data preprocessing complete. Cleaned data saved to Application/data/ ---")
cat("\n")
