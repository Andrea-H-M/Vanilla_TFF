############################################################
## Topic: T-Test for Identifying Functional Signatures    ##
## Author: Olga Andrea Hernandez Miranda, Miranda H       ##
## Date: 02/18/2025                                       ##
## Note: T-Test and Filtering of Functional Signatures    ##
############################################################

# Set working directory
setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutoral3/Enriquesimiento/EnriquesimientoHM/Buscar genes HM/PruebaSignificancia")

# Load required libraries
library(dplyr)

# 1. Read the CSV file containing the data
df <- read.csv("datos_genes_PostPol.csv")

# 2. Check the first few records to ensure data is loaded correctly
head(df)

# 3. Select relevant columns for the 'Post-pol' stage (R1 and R2)
post_pol_data <- df[, c("PostPol_R1", "PostPol_R2")]

# 4. Select only numeric columns for other stages
# Use 'sapply' to filter numeric columns
other_stages_data <- df[, sapply(df, is.numeric) & !names(df) %in% c("PostPol_R1", "PostPol_R2", "genes")]

# 5. Compute the mean expression for each stage
mean_post_pol <- rowMeans(post_pol_data, na.rm = TRUE)  # na.rm = TRUE ignores NA values
mean_other_stages <- rowMeans(other_stages_data, na.rm = TRUE)  # na.rm = TRUE ignores NA values

# 6. Calculate the expression difference between 'Post-pol' and other stages
diff_expression <- mean_post_pol - mean_other_stages

# 7. Add the expression difference to the dataframe
df$diff_expression <- diff_expression

# 8. Perform a statistical analysis comparing 'Post-pol' with other stages

# Conduct the paired t-test for each gene comparing 'Post-pol' vs. other stages
t_test_results <- apply(df[, c("PostPol_R1", "PostPol_R2")], 1, function(x) t.test(x, other_stages_data)$p.value)
t_test_results

# Add t-test results to the dataframe
df$t_test_p_value <- t_test_results

# View the complete dataframe
df

# Initialize a list to store p-values of each comparison
t_test_p_values <- numeric(nrow(df))
df

# Save the dataframe 'df' to a CSV file
write.csv(df, "resultados_prueba_t_PostPol.csv", row.names = FALSE)

# Filter only genes with p-value ≤ 0.00001
filtered_df <- df %>%
  filter(t_test_p_value <= 0.00001)

# View the filtered genes
filtered_df

# Save the filtered dataframe to a CSV file
write.csv(filtered_df, "filtered_genes_PostPol.csv", row.names = FALSE)