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
    log10_ARG_load,
    ARG_div_shan,
    region,
    city,
    sex,
    country,
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
dim(Subset)
```
[1] 60996    14
[1] 60666    14

--> 330 samples removed


### 2. Age sex and and income analyses

Age categogy
```
Subset$sex <- recode(Subset$sex,"female" = "Female", "male"   = "Male")

plot_df <- Subset %>%
  filter(
    !is.na(log10_ARG_load), !is.na(sex), !is.na(age_category)) %>%
  mutate(
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_), Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))) %>%
  filter(!is.na(Income_group))

table_df <- plot_df %>%
  count(Income_group, age_category, sex) %>%
  pivot_wider(
    names_from = sex,
    values_from = n,
    values_fill = 0
  ) %>% mutate(Total = Female + Male)

table_image <- table_df %>%
  gt(groupname_col = "Income_group") %>%
  cols_label(
    age_category = "Age Category",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(title = "Sample Distribution by Age, Sex, and Income Group (HIC vs LMIC)")

table_image
```

Numeric Age
```r
plot_df <- Subset %>% filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(age_years)) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))) %>%filter(!is.na(Income_group))

table_df <- plot_df %>%
  count(Income_group, sex) %>%
  pivot_wider(
    names_from = sex,
    values_from = n,
    values_fill = 0
  ) %>% mutate(Total = Female + Male)

table_image <- table_df %>%
  gt() %>%
  cols_label(
    Income_group = "Income Group",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>% tab_header(title = "Total Sample Counts by Sex and Income Group")

table_image
```

### 3.2 Linear model
Numeric age
```r

HIC_df <- plot_df %>% filter(Income_group == "HIC")
LMIC_df <- plot_df %>% filter(Income_group == "LMIC")

lm_HIC <- lm(log10_ARG_load ~ age_years + sex, data = HIC_df)
summary(lm_HIC)

lm_LMIC <- lm(log10_ARG_load ~ age_years + sex, data = LMIC_df)
summary(lm_LMIC)
```

| Cohort              | Predictor   | Estimate | Std. Error | t value | p-value | Significance |
| ------------------- | ----------- | -------: | ---------: | ------: | ------: | ------------ |
| **HIC (N = 8,257)** | (Intercept) |   2.7303 |    0.00693 |  394.18 |  <2e-16 | ***          |
| HIC                 | Age (years) |  0.00036 |    0.00014 |   2.649 | 0.00809 | **           |
| HIC                 | Sex (Male)  | -0.00263 |    0.00683 |  -0.385 |   0.700 |              |
| **LMIC (N = 986)**  | (Intercept) |   2.7808 |    0.01868 |  148.86 |  <2e-16 | ***          |
| LMIC                | Age (years) |  0.00317 |    0.00037 |   8.483 |  <2e-16 | ***          |
| LMIC                | Sex (Male)  |  0.02930 |    0.01775 |   1.650 |  0.0992 | .            |

| Cohort | Residual SE |      R² | Adjusted R² | F-statistic          | p-value  |
| ------ | ----------- | ------: | ----------: | -------------------- | -------- |
| HIC    | 0.3097      | 0.00086 |     0.00061 | 3.537 (df = 2, 8254) | 0.02914  |
| LMIC   | 0.2787      |  0.0711 |      0.0692 | 37.66 (df = 2, 984)  | <2.2e-16 |


### 3.3 Boxplot
Categorised age
```r
plot_df <- Subset %>% filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(age_category)) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))) %>%filter(!is.na(Income_group))

plot_df <- plot_df %>%
  mutate(
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"
    )),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    sex = recode(sex, "female" = "Female", "male" = "Male"))

ggplot(plot_df, aes(x = age_category, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~Income_group) +  # separate HIC vs LMIC
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by Age Category and Sex",
    x = "Age Category",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  stat_compare_means(
    aes(group = sex), 
    label = "p.signif",    # shows significance symbols
    method = "wilcox.test" # non-parametric comparison
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )
```

### 3.4 LOESS
--> Not statistically best, but nice for visualization.
```
ggplot(plot_df, aes(x = age_years, y = log10_ARG_load, color = sex)) +
  geom_point(alpha = 0.2, size = 1.0) +  # lighter and smaller dots
  geom_smooth(method = "loess", se = TRUE, span = 0.6) +
  facet_wrap(~Income_group) +
  scale_color_npg() +
  labs(
    x = "Age (years)",
    y = "Log10 ARG Load",
    color = "Sex",
    title = "LOESS-smoothed ARG load by age, sex, and income group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA))
```

### 2.6 GAM Female-Male Differences 

```r
df_temp <- Subset |>
  dplyr::select(log10_ARG_load, age_years, sex, World_Bank_Income_Group) |>
  dplyr::mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    age_years = as.numeric(age_years),
    Income_group = dplyr::case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c(
        "Low income",
        "Lower middle income",
        "Upper middle income"
      ) ~ "LMIC",
      TRUE ~ NA_character_
    ),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC"))
  ) |>
  dplyr::filter(
    sex %in% c("Male", "Female"),
    !is.na(Income_group),
    !is.na(age_years))


fit <- gam(log10_ARG_load ~ sex + s(age_years, by = sex, k = 10) + sex*Income_group,
  method = "REML", data= df_temp)

summary(fit)
draw(fit)

```

| Component        | Term                 | Estimate | Std. Error | t / F value | p-value | Significance |
| ---------------- | -------------------- | -------: | ---------: | ----------: | ------: | ------------ |
| **Parametric**   | (Intercept)          |   2.7363 |    0.00482 |      567.42 |  <2e-16 | ***          |
| Parametric       | Sex (Female)         |   0.0103 |    0.00672 |       1.534 |  0.1251 |              |
| Parametric       | Income group (LMIC)  |   0.1959 |    0.01434 |      13.663 |  <2e-16 | ***          |
| Parametric       | Sex (Female) × LMIC  |  -0.0490 |    0.02055 |      -2.383 |  0.0172 | *            |
| **Smooth terms** | s(age_years): Male   |        — |          — |       25.83 |  <2e-16 | ***          |
| Smooth terms     | s(age_years): Female |        — |          — |       19.48 |  <2e-16 | ***          |

| Family   | Link function | Adjusted R² | Deviance explained | REML   | Scale estimate | Sample size (N) |
| -------- | ------------- | ----------- | ------------------ | ------ | -------------- | --------------- |
| Gaussian | Identity      | 6.25%       | 6.41%              | 2095.8 | 0.091347       | 9,244           |

```

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
  sex = c("Female","Male"),
  Income_group = levels(df_temp$Income_group))  # changed from World_Bank_Income_Group)

predictions <- predict(fit, newdata, se.fit = TRUE)
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit

diff_data <- newdata %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = c(pred, se))

diff_data$diff <- diff_data$pred_Male - diff_data$pred_Female

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
  ) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg() +
  theme_minimal()

```

--

## 4. ARG Load of the reproductive age

### 4.1 Create the dataframe
```r
df_income_r_age_ARG <- Subset %>%
  select(
    sex,
    World_Bank_Income_Group,        
    age_years,
    log10_ARG_load
  ) %>%
  mutate(
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c(
        "Low income",
        "Lower middle income",
        "Upper middle income"
      ) ~ "LMIC",
      TRUE ~ NA_character_
    ),
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>%
  filter(
    !is.na(age_years),
    !is.na(sex),
    !is.na(log10_ARG_load),
    !is.na(Income_group),
    (sex == "Female" & age_years >= 15 & age_years <= 49) |
    (sex == "Male" & age_years >= 15 & age_years <= 49)
  )
```

### 4.2 Boxplot
```
ggplot(df_income_r_age_ARG, aes(x = Income_group, y = log10_ARG_load, fill = sex)) +
  geom_jitter(aes(color = sex), width = 0.2, alpha = 0.3, size = 1.5, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    x = "Income Group",
    y = "log10(ARG load)",
    fill = "Sex"
  ) +
  ggtitle("ARG Load by Income Group and Sex\nReproductive ages: Female 15–49, Male 15–49") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
```
### 4.3 Linear model 

```
library(dplyr)

model_HIC <- lm(log10_ARG_load ~ age_years + sex,
  data = df_income_r_age_ARG %>% filter(Income_group == "HIC"))
summary(model_HIC)

model_LMIC <- lm(log10_ARG_load ~ age_years + sex,
  data = df_income_r_age_ARG %>% filter(Income_group == "LMIC"))
summary(model_LMIC)
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

### 2.6 GAM
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
```
### 2.X Model check:
Go through with mahkameh
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

### 2.6 GAM differences
```r
diff <- difference_smooths(
  fit,
  smooth = "s(age_years)",
  n = 200)

draw(diff)
```
### 2.7 Visualize GAM differences
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
```
### Options for visualizations

```r
ggplot(diff_data, aes(age_years, diff)) +
  geom_smooth() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Income_group) +  # changed from World_Bank_Income_Group
  labs(
    y = "Male - Female predicted log10(ARG load)",
    x = "Age (years)",
    title = "Sex Differences in Predicted ARG Load Across Age by Income Group"
  )

ggplot(diff_data, aes(age_years, Income_group, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#b40426",   # red (Female > Male)
    mid = "white",
    high = "#3b4cc0",  # blue (Male > Female)
    midpoint = 0,
    limits = c(-max(abs(diff_data$diff)),
               max(abs(diff_data$diff)))
  ) +
  theme_minimal()

# Blue = Male > Female
# Red = Female > Male
# White = no difference

# 3.Trajectories of predicted means: 

ggplot(newdata, aes(age_years, pred, color = sex, fill = sex)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = pred - 1.96*se,
                  ymax = pred + 1.96*se),
              alpha = 0.2, color = NA) +
  facet_wrap(~Income_group) +
  labs(
    y = "Predicted log10(ARG load)",
    x = "Age (years)",
    title = "Age trajectories of ARG diversity by sex and income group"
  ) +
  theme_minimal()
```
