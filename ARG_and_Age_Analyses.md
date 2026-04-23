## 1. Prepare the data
### 1.1 Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(mgcv)
library(ggpubr)
library(tidyverse)

#dir.create("~/R/library", recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))
.libPaths() #[1] "/users/saankarh/R/library"
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("mia")
# install.packages("gt")
# install.packages("gratia")    

library(gt)
library(gratia)
library(mia)

```
### 1.2 Load TSE
```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

### 1.3 Extract colData
```r
colData_df <- as.data.frame(colData(TSE))

Subset <- colData_df %>%
  select(
    acc,
    log10_ARG_load,
    sex,
    age_years,
    age_category,
    World_Bank_Income_Group,
    country,
    continent,
    parent,
    pregnant,
    bioproject,
    biosample)

Subset[Subset == ""] <- NA
```
### 1.4 Remove samples from Bioproject PRJNA544527
```r
dim(Subset)

remove_ids <- c(
  paste0("SAMN", 11950082:11950286),
  paste0("SAMN", 11950288:11950349),
  paste0("SAMN", 11950351:11950423)
)

Subset <- Subset[!Subset$biosample %in% remove_ids, ]

```
--> 330 samples removed

```

## 2. ARG Load, numeric age, sex and income

### 2.1 Sample distribution
```r
Subset$sex <- recode(Subset$sex, "female" = "Female", "male"   = "Male")

plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load), !is.na(sex), !is.na(age_years))

table_df <- plot_df %>%
  count(sex) %>% pivot_wider(names_from = sex, values_from = n)

table_df <- table_df %>% mutate(Total = Female + Male)

table_image <- table_df %>%
  gt() %>%
  cols_label(
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(title = "Sample Distribution by Sex")

table_image
```


### 2.2 Loess curve of ARG load by sex and numeric age

```r
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$age_years[Subset$age_years == "" | is.na(Subset$age_years)] <- NA

Subset$sex <- recode(Subset$sex, "female" = "Female", "male"   = "Male")

plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(age_years)) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))

ggplot(plot_df, aes(x = age_years, y = log10_ARG_load, color = sex, fill = sex)) +
  geom_point(alpha = 0.08, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7, alpha = 0.2, size = 1.2) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "ARG Load Across Age by Sex",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain")
  )
```

#### 2.4 Linear model
```
n_samples <- Subset %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(age_years)) %>% nrow()

n_samples

lm_full <- lm(log10_ARG_load ~ sex + age_years, data = plot_df)
summary(lm_full)

```

### 2.6 Katan malli

# Differences 

```r
df_temp <- Subset %>%
  select(log10_ARG_load, sex, age_years, World_Bank_Income_Group) %>%
  filter(
    sex %in% c("male", "female"),
    !is.na(age_years), !is.na(World_Bank_Income_Group), !is.na(log10_ARG_load)) %>%
  mutate(
    sex = factor(sex, levels = c("male", "female")),
    age_years = as.numeric(age_years),
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_
    ),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC"))) %>% filter(!is.na(Income_group))

fit <- gam(
  log10_ARG_load ~ sex + s(age_years, by = sex, k = 10) + sex*Income_group,
  method = "REML", data= df_temp)

summary(fit)
draw(fit)

# draw the differences
diff <- difference_smooths(
  fit,
  smooth = "s(age_years)",
  n = 200)

draw(diff)

```
# Model check:
```r
plot(fit$fitted.values, df_temp$log10_ARG_load,
     xlab = "Fitted values",
     ylab = "Observed log10_ARG_load",
     main = "Response vs Fitted Values",
     pch = 16, col = "black")

# Add the diagonal line (y = x)
abline(a = 0, b = 1, col = "red", lwd = 2)

residuals <- df_temp$log10_ARG_load - fit$fitted.values
plot(fit$fitted.values, residuals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h=0, col="red")

# Some exterme values (above 1.5)
# Over 1.0 large resuduals
# 0 to 0.5 normal deviation
```

```r
newdata <- expand.grid(
  age_years = seq(min(df_temp$age_years),
                  max(df_temp$age_years),
                  length = 200),
  sex = c("female","male"),
  Income_group = levels(df_temp$Income_group)  # changed from World_Bank_Income_Group)

predictions <- predict(fit, newdata, se.fit = TRUE)
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit

diff_data <- newdata %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = c(pred, se)
  )


diff_data$diff <- diff_data$pred_male - diff_data$pred_female

# Plot
library(ggplot2)

ggplot(diff_data, aes(age_years, diff)) +
  geom_smooth() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Income_group) +  # changed from World_Bank_Income_Group
  labs(
    y = "Male - Female predicted log10(ARG load)",
    x = "Age (years)",
    title = "Sex Differences in Predicted ARG Load Across Age by Income Group"
  )

```
