
# Load necessary libraries
library(tidyverse)
library(emmeans)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(tidyr)
library(scales)
library(knitr)
library(forcats)
library(lme4)
library(ggeffects)
library(multcomp)
library(performance)
library(car)
library(broom)
library(RColorBrewer)
library(MASS)
library(gt)


# Load your dataset
file_path <- "~/Desktop/Project 4/LeeDataFungalMorphogenesis.csv"
data <- read.csv(file_path)
data <- data %>%
  dplyr::rename(
    Treatment = "Pos..Control..No.Rx...A..A.SA..B..B.SB..combination.of.A.and.B..combination.of.SA.and.SB",
    Replicate = "Replicate",
    ResponseGT = "Response..GT..germtubes..only",
    ResponseGTPH = "Response..GT...PH",
    ResponseGTBuds = "Response..GT...Buds",
    ResponsePHBuds = "Response..PH...Buds",
    ResponseGTMultiple = "Response...1.GTs.multiple.GTs."
  )
# Simplify the treatment names by re-coding them
data$Treatment <- factor(data$Treatment, 
                         levels = c("No Treatment (positive control)", 
                                    "Whole Aloe arborescens Extract",
                                    "compound 'A' isolated from Aloe arborescens",
                                    "compound 'A-SA' from Sigma Aldrich",
                                    "compound 'B' isolated from Aloe arborescens",
                                    "compound 'B-SA' from Sigma Aldrich",
                                    "combination of A and B",
                                    "combination of A-SA and B-SA"),
                         labels = c("Control", 
                                    "WholeExtract", 
                                    "A", 
                                    "A_SA", 
                                    "B", 
                                    "B_SA", 
                                    "A_B", 
                                    "A_SA_B_SA"))
data$Hyphae_Count<- data$ResponseGT + data$ResponseGTPH + data$ResponseGTBuds + data$ResponseGTMultiple 


# Create the Statistical Summary Table
summary_table <- data %>%
  group_by(Treatment) %>%
  summarise(
    Count = n(),
    Mean_Hypha_Count = mean(Hyphae_Count, na.rm = TRUE),
    Median_Hypha_Count = median(Hyphae_Count, na.rm = TRUE),
    SD_Hypha_Count = sd(Hyphae_Count, na.rm = TRUE)
  )


# Convert Mean Hypha Count to Proportion
summary_table <- summary_table %>%
  mutate(Mean_Hypha_Proportion = Mean_Hypha_Count / 200,
         SD_Hypha_Proportion = SD_Hypha_Count / 200)

# Bar Plot of Hypha Count Proportion with Color-Blind-Friendly Palette
bar_plot <- ggplot(summary_table, aes(x = Treatment, y = Mean_Hypha_Proportion, fill = Treatment)) +
  geom_bar(stat = "identity") +  # Use different colors for each treatment
  geom_errorbar(aes(ymin = Mean_Hypha_Proportion - SD_Hypha_Proportion,
                    ymax = Mean_Hypha_Proportion + SD_Hypha_Proportion), width = 0.2) +
  labs(title = "Effect of Treatments on Mean Hypha Formation Proportion",
       x = "Treatment Type",
       y = "Proportion of Hyphae Formed (per 200 cells)") +
  scale_fill_brewer(palette = "Set2") +               # Use a color-blind-friendly palette
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),           # Center the title
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.text = element_text(size = 10),              # Adjust axis label size
    legend.position = "none",                         # Remove the legend
    panel.grid = element_blank()                     # Remove the grid backgrou  # Add axis border
  )


data$offset_log200 <- log(200)


# Fit the quasi-Poisson model
data$offset_log200 <- log(200)
model <- glm(Hyphae_Count ~ Treatment + offset(offset_log200), 
             family = quasipoisson(link = "log"), data = data)




# Fit the quasi-Poisson model
data$offset_log200 <- log(200)
model <- glm(Hyphae_Count ~ Treatment + offset(offset_log200), 
             family = quasipoisson(link = "log"), data = data)

# Extract important summary statistics with confidence intervals
summary_table <- tidy(model, conf.int = TRUE) %>%
  mutate(
    Estimate = exp(estimate),                           # Exponentiated effect size
    `Lower 95% CI` = exp(conf.low),                     # Exponentiate lower CI
    `Upper 95% CI` = exp(conf.high),                    # Exponentiate upper CI
    Significance = case_when(                           # Significance codes
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(term, Estimate, `Lower 95% CI`, `Upper 95% CI`, std.error, statistic, p.value, Significance) %>%
  dplyr::rename(
    "Treatment" = term,
    "Std. Error" = std.error,
    "t value" = statistic,
    "P-value" = p.value
  )



# Calculate the dispersion statistic
dispersion_stat <- summary(model)$deviance / summary(model)$df.residual


# Extract residuals
model_data <- data.frame(
  Fitted = fitted(model),
  Pearson_Residuals = residuals(model, type = "pearson"),
  Deviance_Residuals = residuals(model, type = "deviance")
)

# Plot Pearson residuals
g<-ggplot(model_data, aes(x = Fitted, y = Pearson_Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Pearson Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Pearson Residuals") +
  theme_minimal()

# Plot Deviance residuals
p<-ggplot(model_data, aes(x = Fitted, y = Deviance_Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Deviance Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Deviance Residuals") +
  theme_minimal()
# Goodness-of-fit test
gof.pvalue = 1 - pchisq(model$deviance, model$df.residual)



# Obtain estimated marginal means and set up contrasts for each treatment vs control
emm <- emmeans(model, ~ Treatment)
contrast1 <- contrast(emm, method = "trt.vs.ctrl", ref = "Control")

# Get summary results for Hypothesis 1
contrast1_summary <- summary(contrast1, infer = TRUE)

# Convert summary to a data frame for formatting with gt
contrast1_df <- as.data.frame(contrast1_summary)

# Create a formatted table with gt (adjust columns based on actual names if needed)
gt_table <- contrast1_df %>%
  gt() %>%
  cols_label(
    contrast = "Comparison",
    estimate = "Estimate",
    SE = "Std. Error",
    asymp.LCL = "95% CI Lower",
    asymp.UCL = "95% CI Upper",
    z.ratio = "Z-ratio",
    p.value = "P-value"
  ) %>%
  fmt_number(
    columns = vars(estimate, SE, asymp.LCL, asymp.UCL, z.ratio),
    decimals = 2
  ) %>%
  fmt_number(
    columns = vars(p.value),
    decimals = 4,
    use_seps = FALSE
  ) %>%
  tab_header(
    title = "Hypothesis 1: Each Treatment vs. Control"
  )



# Fit the quasi-Poisson model with an offset
data$offset_log200 <- log(200)
model <- glm(Hyphae_Count ~ Treatment + offset(offset_log200), 
             family = quasipoisson(link = "log"), data = data)

# Obtain estimated marginal means for each treatment
emm <- emmeans(model, ~ Treatment)

# Define custom contrasts for Hypothesis 2 as a list
contrast2 <- contrast(emm, method = list(
  "A vs A_SA" = c("Control" = 0, "WholeExtract" = 0, "A" = 1, "A_SA" = -1, "B" = 0, "B_SA" = 0, "A_B" = 0, "A_SA_B_SA" = 0),
  "B vs B_SA" = c("Control" = 0, "WholeExtract" = 0, "A" = 0, "A_SA" = 0, "B" = 1, "B_SA" = -1, "A_B" = 0, "A_SA_B_SA" = 0),
  "A_B vs A_SA_B_SA" = c("Control" = 0, "WholeExtract" = 0, "A" = 0, "A_SA" = 0, "B" = 0, "B_SA" = 0, "A_B" = 1, "A_SA_B_SA" = -1)
))

# Get summary results for Hypothesis 2 with Bonferroni adjustment
contrast2_summary <- summary(contrast2, infer = TRUE, adjust = "bonferroni")

# Convert summary to a data frame for formatting with gt
contrast2_df <- as.data.frame(contrast2_summary)

# Create a nicely formatted table with gt
gt_table <- contrast2_df %>%
  gt() %>%
  cols_label(
    contrast = "Comparison",
    estimate = "Estimate",
    SE = "Std. Error",
    asymp.LCL = "95% CI Lower",
    asymp.UCL = "95% CI Upper",
    z.ratio = "Z-ratio",
    p.value = "P-value (Bonferroni)"
  ) %>%
  fmt_number(
    columns = vars(estimate, SE, asymp.LCL, asymp.UCL, z.ratio),
    decimals = 2
  ) %>%
  fmt_number(
    columns = vars(p.value),
    decimals = 4,
    use_seps = FALSE
  ) %>%
  tab_header(
    title = "Hypothesis 2: Specific Group Comparisons with Bonferroni Adjustment"
  )



# Fit the quasi-Poisson model with an offset
data$offset_log200 <- log(200)
model <- glm(Hyphae_Count ~ Treatment + offset(offset_log200), 
             family = quasipoisson(link = "log"), data = data)

# Obtain estimated marginal means for each treatment
emm <- emmeans(model, ~ Treatment)

# Set up contrasts comparing each treatment to WholeExtract as the base
contrast3 <- contrast(emm, method = "trt.vs.ctrl", ref = "WholeExtract")

# Get summary results for Hypothesis 3
contrast3_summary <- summary(contrast3, infer = TRUE)

# Convert summary to a data frame for formatting with gt
contrast3_df <- as.data.frame(contrast3_summary)

# Create a nicely formatted table with gt
gt_table <- contrast3_df %>%
  gt() %>%
  cols_label(
    contrast = "Comparison",
    estimate = "Estimate",
    SE = "Std. Error",
    asymp.LCL = "95% CI Lower",
    asymp.UCL = "95% CI Upper",
    z.ratio = "Z-ratio",
    p.value = "P-value"
  ) %>%
  fmt_number(
    columns = vars(estimate, SE, asymp.LCL, asymp.UCL, z.ratio),
    decimals = 2
  ) %>%
  fmt_number(
    columns = vars(p.value),
    decimals = 4,
    use_seps = FALSE
  ) %>%
  tab_header(
    title = "Hypothesis 3: Each Treatment vs. WholeExtract"
  )
