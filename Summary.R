library(mgcv)
library(readxl)
library(ggeffects)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(nlme)
library(visreg)
library(gmodels)
library(metafor)
source("MMK.R")


# ---------------------------------------------------
# NOTE
# ---------------------------------------------------
#This summary is based on the “time series from Denmark” dataset (the third script, using time_series.xlsx and crosstable_DK.xlsx). The same workflow applies for the other two time series datasets—Netherlands (Netherland.xlsx, crosstable_NL.xlsx) and Sweden (Sweden.xlsx, crosstable_SW.xlsx)—with a few minor modifications:
#Data file names: Update read_excel() to load the appropriate country-specific files.
#Column selection: Adjust the column ranges in select() (e.g., 2:141 for Netherlands, 2:71 for Sweden, 2:130 for Denmark) to cover the correct number of taxon columns.
#Plotting year and richness ranges:
#Modify scale_x_continuous(limits = ..., breaks = ...) to fit the actual year range in each dataset.
#Adjust scale_y_continuous(limits = ...) if species richness values differ between countries.
#All other code and processing steps remain the same; only these dataset-specific details need to be updated for each time series.

# ---------------------------------------------------
# Import and clean data
# ---------------------------------------------------
mydata <- as.data.frame(read_excel("dataset.xlsx"))
mydata <- mydata[!grepl("/", mydata$taxon), ]
mydata <- mydata[!grepl("sp.", mydata$taxon), ]
mydata <- mydata[!grepl("spp.", mydata$taxon), ]

Annual <- mydata %>%
  group_by(year) %>%
  summarize(c(Species = length(taxon), sites = length(site_id)))
Annual

# ---------------------------------------------------
# Fit without missing values
# ---------------------------------------------------
set.seed(123)
mydata2 <- as.data.frame(read_excel("crosstable_DK.xlsx"))
mydata22 <- mydata2 %>% select(2:130)
mydata22$Richness <- rowSums(mydata22 > 0)
mydata23 <- mydata2 %>% select(2:130)
mydata23$Total_Abundance <- apply(mydata23, 1, function(x) sum(x[which(x > 0)]))
df_no_gap <- cbind(mydata2$year, mydata2$iteration, mydata22$Richness, mydata23$Total_Abundance)
colnames(df_no_gap) <- c("year", "richness", "tot_abundance")
df_no_gap <- as.data.frame(df_no_gap)

df_no_gap <- filter(df_no_gap, !is.na(richness))
Richness_0missing <- df_no_gap %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(mean_richness = mean(richness), SD_richness = sd(richness))
colnames(Richness_0missing) <- c("year", "mean_richness", "sd_rich")

plot0miss <- ggplot(Richness_0missing, aes(x = year, y = mean_richness)) +
  geom_point(color = "blue", size = 2) +
  scale_x_continuous(limits = c(1990, 2020), breaks = seq(1990, 2020, 5), labels = seq(1990, 2020, 5)) +
  scale_y_continuous(limits = c(20, 50), breaks = seq(20, 50, 5), labels = seq(20, 50, 5)) +
  geom_smooth() +
  theme_classic2() +
  theme_cleveland()
plot0miss

gam_0miss <- gam(mean_richness ~ s(year), data = Richness_0missing, method = "REML")
fit <- as.data.frame(gam_0miss$fitted.values)
AIC(gam_0miss)
colnames(fit)[1] <- "annual_fit"
years <- data.frame(year = seq(1992, 2018))
original_full_fit <- cbind(years, fit)

rmse_orig_gam <- sqrt(mean(residuals(gam_0miss)^2))

lm_0miss <- lm(mean_richness ~ year, data = Richness_0missing)
fit <- as.data.frame(lm_0miss$fitted.values)
colnames(fit)[1] <- "annual_fit"
original_full_fit_lm <- cbind(years, fit)
rmse_orig_lm <- sqrt(mean(residuals(lm_0miss)^2))

# ---------------------------------------------------
# 1 gap (deleting 1 row randomly, 100 iterations)
# ---------------------------------------------------
set.seed(123)
mydata2 <- as.data.frame(read_excel("crosstable_DK.xlsx"))
final_df_1_missing <- data.frame()

for (i in 1:100) {
  random_row <- sample(2:(nrow(mydata2) - 1), 1)
  first_col_value <- mydata2[random_row, 1]
  mydata3 <- mydata2[-random_row, ]
  mydata3 <- rbind(mydata3, c(first_col_value, rep(NA, ncol(mydata3) - 1)))
  mydata3 <- mydata3[order(mydata3[, 1]), ]
  mydata3$iteration <- i
  final_df_1_missing <- rbind(final_df_1_missing, mydata3)
}

final_df_1_missing2 <- final_df_1_missing %>% select(2:130)
final_df_1_missing2$Richness <- rowSums(final_df_1_missing2 > 0)
final_df_1_missing3 <- final_df_1_missing %>% select(2:130)
final_df_1_missing3$Total_Abundance <- apply(final_df_1_missing3, 1, function(x) sum(x[which(x > 0)]))
final_df_1_gap <- cbind(final_df_1_missing$year, final_df_1_missing$iteration, final_df_1_missing2$Richness, final_df_1_missing3$Total_Abundance)
colnames(final_df_1_gap) <- c("year", "iteration", "richness", "tot_abundance")
final_df_1_gap <- as.data.frame(final_df_1_gap)

final_df_1_gap <- filter(final_df_1_gap, !is.na(richness))
Richness_1missing <- final_df_1_gap %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(mean_richness = mean(richness), SD_richness = sd(richness))
colnames(Richness_1missing) <- c("year", "mean_richness", "sd_rich")

# Plot and GAM/LM fits for 1 gap
plot1miss <- ggplot(Richness_1missing, aes(x = year, y = mean_richness)) +
  geom_point(color = "blue", size = 2) +
  scale_x_continuous(limits = c(1990, 2020), breaks = seq(1990, 2020, 5), labels = seq(1990, 2020, 5)) +
  scale_y_continuous(limits = c(20, 50), breaks = seq(0, 50, 5), labels = seq(0, 50, 5)) +
  geom_smooth() +
  theme_classic2() +
  theme_cleveland()
plot1miss

gam_1miss <- gam(mean_richness ~ s(year), data = Richness_1missing, method = "REML")
fit <- as.data.frame(gam_1miss$fitted.values)
colnames(fit)[1] <- "annual_fit"
years <- data.frame(year = seq(1992, 2018))
fit1miss <- cbind(years, fit)

# Mann-Kendall, GAM and LM loop for 1 gap
MK_1_miss <- as.data.frame(matrix(nrow = 0, ncol = 11))
colnames(MK_1_miss) <- c("Corrected Zc", "new P-value", "N/N*", "Original Z", "old P.value", "Tau", "Sen s slope", "old.variance", "new.variance", "S statistic", "n")

iteration <- unique(final_df_1_gap$iteration)
for (n in iteration) {
  subset_df <- final_df_1_gap[final_df_1_gap$iteration == n, ]
  MK <- My.mmkh(subset_df$richness)
  MK_1_miss <- rbind(MK_1_miss, MK)
  cat(".")
}

AIC_1_miss_gam <- data.frame(AIC = numeric(0))
DEV_1_miss_gam <- data.frame(DEV = numeric(0))
RMSE_1_miss_gam <- data.frame(DEV = numeric(0))
loop_1_missing <- data.frame(fit = numeric(0))

for (n in iteration) {
  subset_df <- final_df_1_gap[final_df_1_gap$iteration == n, ]
  gam_1miss2 <- gam(richness ~ s(year), data = subset_df)
  fit <- gam_1miss2$fitted.values
  AIC_1_miss_gam <- rbind(AIC_1_miss_gam, data.frame(AIC = AIC(gam_1miss2)))
  rmse <- sqrt(mean(residuals(gam_1miss2)^2))
  RMSE_1_miss_gam <- rbind(RMSE_1_miss_gam, data.frame(DEV = rmse))
  dev_percent <- ((gam_1miss2$null.deviance - gam_1miss2$deviance) / gam_1miss2$null.deviance) * 100
  DEV_1_miss_gam <- rbind(DEV_1_miss_gam, data.frame(DEV = dev_percent))
  loop_1_missing <- rbind(loop_1_missing, data.frame(fit = fit))
  cat(".")
}
fit1miss_gam <- cbind(final_df_1_gap$year, final_df_1_gap$iteration, loop_1_missing)
colnames(fit1miss_gam) <- c("year", "iteration", "fit")

AIC_1_miss_lm <- data.frame(AIC = numeric(0))
DEV_1_miss_lm <- data.frame(DEV = numeric(0))
RMSE_1_miss_lm <- data.frame(DEV = numeric(0))
loop_1_missing_lm <- data.frame(fit = numeric(0))

for (n in iteration) {
  subset_df <- final_df_1_gap[final_df_1_gap$iteration == n, ]
  lm_1miss2 <- lm(richness ~ year, data = subset_df)
  fit <- lm_1miss2$fitted.values
  AIC_1_miss_lm <- rbind(AIC_1_miss_lm, data.frame(AIC = AIC(lm_1miss2)))
  rmse <- sqrt(mean(residuals(lm_1miss2)^2))
  RMSE_1_miss_lm <- rbind(RMSE_1_miss_lm, data.frame(DEV = rmse))
  dev <- deviance(lm_1miss2) / lm_1miss2$df.residual
  DEV_1_miss_lm <- rbind(DEV_1_miss_lm, data.frame(DEV = dev))
  loop_1_missing_lm <- rbind(loop_1_missing_lm, data.frame(fit = fit))
  cat(".")
}
fit1miss_lm <- cbind(final_df_1_gap$year, final_df_1_gap$iteration, loop_1_missing_lm)
colnames(fit1miss_lm) <- c("year", "iteration", "fit")

# Plots for 1 gap (raw, LM, GAM) here...

# ---------------------------------------------------
# 2 gaps, 5 gaps, 10 gaps
# ---------------------------------------------------
# For "2 gaps", "5 gaps", and "10 gaps", repeat the exact same process as for "1 gap" above,
# but instead of deleting one random row, delete 2, 5, or 10 random rows for each iteration.
# For each scenario:
#   - Create the "final_df_x_missing" by deleting x random rows per iteration (x = 2, 5, 10)
#   - Calculate "Richness_xmissing", run the same MK, GAM, and LM loops as above
#   - Gather the fits (fitxmiss_gam, fitxmiss_lm)
#   - Plot the results for each scenario
#   - All code structure and plot types remain the same, just replace "1" with the correct gap number

# ---------------------------------------------------
# Final combined plots
# ---------------------------------------------------
# At the end, you generate:
#   - A grid plot comparing LM and GAM fits for 1, 2, 5, 10 gaps side-by-side
#   - Overlay plots comparing the fits for all gap scenarios together
#   - Raw data overlays for all scenarios
#   - Export plots as SVG with cowplot::plot_grid()

library(cowplot)
plot1 <- plot_grid(
  plot1miss_lm, plot1miss_gam,
  plot2miss_lm, plot2miss_gam,
  plot5miss_lm, plot5miss_gam,
  plot10miss_lm, plot10miss_gam,
  ncol = 2
)
svg("plot_ALL_missing_DK_sep.svg")
plot1
dev.off()

plot2 <- plot_grid(
  plot_miss_overall_LM,
  plot_miss_overall_GAM,
  plot_miss_overall_raw,
  ncol = 1
)
svg("plot_ALL_missing_DK.svg")
plot2
dev.off()

# ---------------------------------------------------
# Exporting average and SD of AICs for each gap and model
# ---------------------------------------------------
mean(AIC_1_miss_lm$AIC); sd(AIC_1_miss_lm$AIC)
mean(AIC_2_miss_lm$AIC); sd(AIC_2_miss_lm$AIC)
mean(AIC_5_miss_lm$AIC); sd(AIC_5_miss_lm$AIC)
mean(AIC_10_miss_lm$AIC); sd(AIC_10_miss_lm$AIC)
mean(AIC_1_miss_gam$AIC); sd(AIC_1_miss_gam$AIC)
mean(AIC_2_miss_gam$AIC); sd(AIC_2_miss_gam$AIC)
mean(AIC_5_miss_gam$AIC); sd(AIC_5_miss_gam$AIC)
mean(AIC_10_miss_gam$AIC); sd(AIC_10_miss_gam$AIC)

# ---------------------------------------------------
# Data transformation: crosstable to long format for 1, 2, 5, 10 gaps
# ---------------------------------------------------
library(tidyr)
final_df_1_missing_list  <- pivot_longer(final_df_1_missing,  cols = !c("year", "iteration"),
                                         names_to = "taxon", values_to = "abundance")
final_df_2_missing_list  <- pivot_longer(final_df_2_missing,  cols = !c("year", "iteration"),
                                         names_to = "taxon", values_to = "abundance")
final_df_5_missing_list  <- pivot_longer(final_df_5_missing,  cols = !c("year", "iteration"),
                                         names_to = "taxon", values_to = "abundance")
final_df_10_missing_list <- pivot_longer(final_df_10_missing, cols = !c("year", "iteration"),
                                         names_to = "taxon", values_to = "abundance")

# ---------------------------------------------------
# MICE multiple imputation loop (1 gap, shown in full)
# ---------------------------------------------------
library(mice)
mice_1_filled <- data.frame()
iteration <- unique(final_df_1_missing_list$iteration)

for(n in iteration){
  subset_df <- final_df_1_missing_list[final_df_1_missing_list$iteration == n, ]
  species <- unique(subset_df$taxon)
  for(j in 1:length(species)){
    subset_species_df <- subset_df %>% filter(taxon == species[j])
    if(sum(!is.na(subset_species_df$abundance)) < 2) next
    new_results <- complete(mice(subset_species_df, m = 5, maxit = 50, method = "pmm", seed = 123))
    mice_1_filled <- rbind(mice_1_filled, new_results)
    cat("Imputation for iteration", n, "species", species[j], "\n")
  }
}

# ---------------------------------------------------
# For 2, 5, 10 gaps: repeat same MICE loop structure as above, 
# but use final_df_2_missing_list, final_df_5_missing_list, final_df_10_missing_list
# and output to mice_2_filled, mice_5_filled, mice_10_filled.
# ---------------------------------------------------
# The process for 2, 5, 10 gaps is identical: 
#  - For each iteration, for each species, run mice imputation as above
#  - Store results in mice_x_filled dataframes

# ---------------------------------------------------
# Cleaning NAs after MICE for all filled dataframes
# ---------------------------------------------------
mice_1_filled[is.na(mice_1_filled)]     <- 0
mice_2_filled[is.na(mice_2_filled)]     <- 0
mice_5_filled[is.na(mice_5_filled)]     <- 0
mice_10_filled[is.na(mice_10_filled)]   <- 0

# ---------------------------------------------------
# Convert imputed long data back to wide format and calculate richness/abundance
# (Detailed for 1 gap below)
# ---------------------------------------------------
mice_1_filled$iteration <- as.character(mice_1_filled$iteration)
test1 <- pivot_wider(mice_1_filled, id_cols = c("year", "iteration"),
                     names_from = "taxon", values_from = "abundance")
test12 <- test1 %>% select(3:131)
test12$Richness <- rowSums(test12 > 0)
test13 <- test1 %>% select(3:131)
test13$Total_Abundance <- apply(test13, 1, function(x) sum(x[x > 0]))
final_df_1_filled <- cbind(test1$year, test1$iteration, test12$Richness, test13$Total_Abundance)
colnames(final_df_1_filled) <- c("year", "iteration", "richness", "tot_abundance")
final_df_1_filled <- as.data.frame(final_df_1_filled)

# ---------------------------------------------------
# Repeat the wide-formatting and richness/abundance calculation 
# for mice_2_filled, mice_5_filled, mice_10_filled.
# Output to: final_df_2_filled, final_df_5_filled, final_df_10_filled.
# ---------------------------------------------------

# ---------------------------------------------------
# Calculate annual species richness + SD for each filled dataset
# ---------------------------------------------------
final_df_1_filled <- filter(final_df_1_filled, !is.na(richness))
final_df_1_filled$year <- as.numeric(final_df_1_filled$year)
final_df_1_filled$iteration <- as.numeric(final_df_1_filled$iteration)
final_df_1_filled$richness <- as.numeric(final_df_1_filled$richness)

# Repeat above for 2, 5, 10 gaps if needed

Richness_filled <- final_df_1_filled %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(mean_richness = mean(richness), SD_richness = sd(richness))
colnames(Richness_filled) <- c("year", "mean_richness", "sd_rich")
Richness_filled <- as.data.frame(Richness_filled)

# ---------------------------------------------------
# Plot results after filling
# ---------------------------------------------------
plot10miss <- ggplot(Richness_filled, aes(x = year, y = mean_richness)) +
  geom_point(color = "blue", size = 2) +
  scale_x_continuous(limits = c(1990, 2020), breaks = seq(1990, 2020, 5), labels = seq(1990, 2020, 5)) +
  scale_y_continuous(limits = c(20, 50), breaks = seq(0, 50, 5), labels = seq(0, 50, 5)) +
  geom_smooth() +
  theme_classic2() +
  theme_cleveland()
plot10miss

gam_filled <- gam(mean_richness ~ s(year), data = Richness_filled, method = "REML")
fit <- as.data.frame(gam_filled$fitted.values)
colnames(fit)[1] <- "annual_fit"
years <- data.frame(year = seq(1992, 2018))
fit1_filled <- cbind(years, fit)

# FILLED MODELS

# --- 1 filled/gap: FULL CODE ---

# MK trend for 1 filled
MK_1_fill <- as.data.frame(matrix(nrow=0,ncol=11))
columns = c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_1_fill) <- columns

iteration <- unique(final_df_1_filled$iteration)
for(n in iteration){
  subset_df <- final_df_1_filled[final_df_1_filled$iteration == n,]
  MK <- My.mmkh(subset_df$richness)
  MK_1_fill <- rbind(MK_1_fill, MK)
  names(MK_1_fill) <- columns
  cat(".")
}

# GAM model for 1 filled
AIC_1_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "AIC"; colnames(AIC_1_filled_gam) <- columns
DEV_1_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "DEV"; colnames(DEV_1_filled_gam) <- columns
RMSE_1_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "DEV"; colnames(RMSE_1_filled_gam) <- columns
loop_1_filled <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "fit"; colnames(loop_1_filled) <- columns

iteration <- unique(final_df_1_filled$iteration)
for(n in iteration){
  subset_df <- final_df_1_filled[final_df_1_filled$iteration == n,]
  gam_1_filled <- gam(richness ~ s(year), data = subset_df)
  fit <- gam_1_filled$fitted.values
  AIC <- as.data.frame(AIC(gam_1_filled)); colnames(AIC)[1] <- "AIC"; AIC_1_filled_gam <- rbind(AIC_1_filled_gam, AIC)
  rmse <- sqrt(mean(residuals(gam_1_filled)^2)); RMSE_1_filled_gam <- rbind(RMSE_1_filled_gam, rmse)
  DEV <- ((gam_1_filled$null.deviance - gam_1_filled$deviance)/gam_1_filled$null.deviance)*100
  DEV <- as.data.frame(DEV); colnames(DEV)[1] <- "DEV"; DEV_1_filled_gam <- rbind(DEV_1_filled_gam, DEV)
  fit <- as.data.frame(fit); colnames(fit)[1] <- "fit"; loop_1_filled <- rbind(loop_1_filled, fit)
  cat(".")
}
fit1_filled_gam <- cbind(final_df_1_filled$year, final_df_1_filled$iteration, loop_1_filled)
colnames(fit1_filled_gam) <- c("year", "iteration", "fit")

# Linear model for 1 filled
AIC_1_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "AIC"; colnames(AIC_1_filled_lm) <- columns
DEV_1_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "DEV"; colnames(DEV_1_filled_lm) <- columns
RMSE_1_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "DEV"; colnames(RMSE_1_filled_lm) <- columns
loop_1_filled <- as.data.frame(matrix(nrow=0,ncol=1)); columns = "fit"; colnames(loop_1_filled) <- columns

iteration <- unique(final_df_1_filled$iteration)
for(n in iteration){
  subset_df <- final_df_1_filled[final_df_1_filled$iteration == n,]
  lm_1_filled <- lm(richness ~ year, data = subset_df)
  fit <- lm_1_filled$fitted.values
  AIC <- as.data.frame(AIC(lm_1_filled)); colnames(AIC)[1] <- "AIC"; AIC_1_filled_lm <- rbind(AIC_1_filled_lm, AIC)
  rmse <- sqrt(mean(residuals(lm_1_filled)^2)); RMSE_1_filled_lm <- rbind(RMSE_1_filled_lm, rmse)
  DEV <- deviance(lm_1_filled)/lm_1_filled$df.residual; DEV <- as.data.frame(DEV); colnames(DEV)[1] <- "DEV"; DEV_1_filled_lm <- rbind(DEV_1_filled_lm, DEV)
  fit <- as.data.frame(fit); colnames(fit)[1] <- "fit"; loop_1_filled <- rbind(loop_1_filled, fit)
  cat(".")
}
fit1_filled_lm <- cbind(final_df_1_filled$year, final_df_1_filled$iteration, loop_1_filled)
colnames(fit1_filled_lm) <- c("year", "iteration", "fit")

# Plots for 1 filled
plot1fill_lm <- ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(30, 40,5), labels = seq(30, 40,5))+
  geom_smooth(data=fit1_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  theme_classic2()+ theme_cleveland()

plot1fill_gam <- ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit1_filled_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "gam", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "gam",se=FALSE,colour="red")+
  theme_classic2()+ theme_cleveland()

plot1fill_raw <- ggplot()+
  geom_smooth(data=final_df_1_gap, aes(x=year, y=richness, group=iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_point(data=final_df_1_gap,aes(x=year, y=richness, group=iteration),color = "blue", size = 3)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+ theme_cleveland()

# --- 2, 5, and 10 filled/gaps: SHORT INSTRUCTIONS ---

# For 2 filled/gap, copy the code for "1 filled", 
# but change all instances of "1" to "2" in variable names and in the dataframe used (e.g., use final_df_2_filled).
# Do the same for MK_2_fill, fit2_filled_gam, fit2_filled_lm, plot2fill_lm, plot2fill_gam, etc.

# For 5 filled/gap, repeat as above but change "1" to "5" and use the relevant dataframe final_df_5_filled, 
# and so on for all related variables and plots.

# For 10 filled/gap, again, copy and modify as above, using final_df_10_filled and so on. 
# Each section follows exactly the same structure as for 1 filled/gap, with only the numeric identifier and data source changed.

# Exporting average and standard deviation of AICs for 1 filled and 1 missing

# Linear model (LM) for 1 filled
mean(AIC_1_filled_lm$AIC)
sd(AIC_1_filled_lm$AIC)

# Generalized additive model (GAM) for 1 filled
mean(AIC_1_filled_gam$AIC)
sd(AIC_1_filled_gam$AIC)

# Repeat the same for 2, 5, and 10 filled or missing datasets:
# - Copy the above lines, but change the "1" in variable names to "2", "5", and "10" respectively.
# - Use AIC_2_filled_lm, AIC_5_filled_lm, etc., for LM,
#   and AIC_2_filled_gam, AIC_5_filled_gam, etc., for GAM.
# - Do this for both "filled" and "miss" as needed.

# ------------------------------------------------------

# Averaging the richness values before AIC (example for 1 gap/missing and 1 filled)

# LM: 1 missing (gaps)
Annual <- as.data.frame(final_df_1_gap %>%
                          group_by(year) %>%
                          summarize(mean = mean(richness)))
avg1_miss <- lm(mean ~ year, data = Annual)
AIC(avg1_miss)

# LM: 1 filled
Annual <- as.data.frame(final_df_1_filled %>%
                          group_by(year) %>%
                          summarize(mean = mean(richness)))
avg1_fill <- lm(mean ~ year, data = Annual)
AIC(avg1_fill)

# GAM: 1 missing (gaps)
Annual <- as.data.frame(final_df_1_gap %>%
                          group_by(year) %>%
                          summarize(mean = mean(richness)))
avg1_miss_gam <- gam(mean ~ s(year), data = Annual)
AIC(avg1_miss_gam)

# GAM: 1 filled
Annual <- as.data.frame(final_df_1_filled %>%
                          group_by(year) %>%
                          summarize(mean = mean(richness)))
avg1_fill_gam <- gam(mean ~ s(year), data = Annual)
AIC(avg1_fill_gam)

# To do the same for 2, 5, and 10 gaps/filled:
# - Change all "1" to "2", "5", or "10" in the variable and dataframe names.
# - For example, use final_df_2_gap, final_df_2_filled, etc.
# - Resulting models are avg2_miss, avg5_miss, avg10_miss, etc.

# To create summary tables:
# Combine AICs using rbind, label columns and rows accordingly:
AICs_missing_lm_DK <- rbind(AIC(avg1_miss), AIC(avg2_miss), AIC(avg5_miss), AIC(avg10_miss))
colnames(AICs_missing_lm_DK) <- c("AIC_lm_missing")
rownames(AICs_missing_lm_DK) <- c("1","2","5","10")
# Repeat this rbind approach for filled and for GAM models.

# ------------------------------------------------------

# MK Stuff - Renaming columns for 1 gap/fill
colnames(MK_1_miss)[10] <- "S_stat"
colnames(MK_1_fill)[10] <- "S_stat"
# Repeat for MK_2_miss, MK_5_miss, MK_10_miss, MK_2_fill, etc.
# Just swap the numeric identifier accordingly.

# ------------------------------------------------------

# Meta regression and meta-analysis plots for 1 gap/fill

orig <- My.mmkh(Richness_0missing$mean_richness)
orig <- as.data.frame(t(orig))
colnames(orig)[10] <- "S_stat"
MK_orig <- rma(S_stat, old.variance, data = orig)

MK_trend_1_miss <- rma(S_stat, old.variance, data = MK_1_miss)
MK_trend_1_fill <- rma(S_stat, old.variance, data = MK_1_fill)

# To do for 2, 5, 10:
# - Copy the code for MK_trend_1_miss and MK_trend_1_fill, 
# - Change all "1" to "2", "5", or "10" in variable and dataframe names.

# Combine meta-regression results into tables:
t0 <- cbind(MK_orig$b, MK_orig$se, MK_orig$ci.lb, MK_orig$ci.ub)
t1 <- cbind(MK_trend_1_miss$b, MK_trend_1_miss$se, MK_trend_1_miss$ci.lb, MK_trend_1_miss$ci.ub)
q1 <- cbind(MK_trend_1_fill$b, MK_trend_1_fill$se, MK_trend_1_fill$ci.lb, MK_trend_1_fill$ci.ub)
# For 2, 5, 10 repeat as above.

pl_miss_DK <- rbind(t0, t1) # add t2, t3, t4 for 2, 5, 10
rownames(pl_miss_DK) <- c("0","1","2","5","10")
pl_miss_DK <- cbind(rownames(pl_miss_DK), data.frame(pl_miss_DK, row.names=NULL))
colnames(pl_miss_DK) <- c("rep","estimate","se","ci.lb","ci.ub")

pl_fill_DK <- rbind(t0, q1) # add q2, q3, q4 for 2, 5, 10
rownames(pl_fill_DK) <- c("0","1","2","5","10")
pl_fill_DK <- cbind(rownames(pl_fill_DK), data.frame(pl_fill_DK, row.names=NULL))
colnames(pl_fill_DK) <- c("rep","estimate","se","ci.lb","ci.ub")

# Plotting the meta-regression results
MK_plot_miss <- ggplot(pl_miss_DK, aes(x=rep, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub, width=0.5)) +
  theme_bw() + theme_cleveland()

MK_plot_fill <- ggplot(pl_fill_DK, aes(x=rep, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub, width=0.5)) +
  theme_bw() + theme_cleveland()

# ------------------------------------------------------

# Individual model plots (example for 1 gap/fill)
library(sjPlot)
p <- plot_model(MK_trend_1_miss, type = "est")
pMK_miss_1 <- p + theme_sjplot()

p <- plot_model(MK_trend_1_fill, type = "est")
pMK_fill_1 <- p + theme_sjplot()

# For 2, 5, 10: Repeat as above, changing the identifier.

# Combine plots as needed with plot_grid (from cowplot), as above.

# ------------------------------------------------------

# RMSE calculations (example for 1 gap/fill)

colnames(RMSE_1_miss_gam)[1] <- "value"
colnames(RMSE_1_filled_gam)[1] <- "value"
colnames(RMSE_1_miss_lm)[1] <- "value"
colnames(RMSE_1_filled_lm)[1] <- "value"

mean_1_miss_gam <- mean(RMSE_1_miss_gam$value)
se_1_miss_gam <- sd(RMSE_1_miss_gam$value)
mean_1_fill_gam <- mean(RMSE_1_filled_gam$value)
se_1_fill_gam <- sd(RMSE_1_filled_gam$value)
mean_1_miss_lm <- mean(RMSE_1_miss_lm$value)
se_1_miss_lm <- sd(RMSE_1_miss_lm$value)
mean_1_fill_lm <- mean(RMSE_1_filled_lm$value)
se_1_fill_lm <- sd(RMSE_1_filled_lm$value)

# For 2, 5, 10, repeat the above replacing "1" with "2", "5", and "10" as needed.

# Combine means/SEs for summary tables as above, and plot using ggplot.

# ------------------------------------------------------

# In summary:
# For all code sections above, after writing for 1 gap/fill, 
# repeat the same logic for 2, 5, and 10 by swapping the identifiers in variable and dataframe names.
# The code structure remains exactly the same.

# This approach keeps your script DRY and clear!


