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
```
--> 330 samples removed

--

## 2. Age sex and and income analyses

### 2.1 Sample distriburions

Age categogy
```r
Subset$sex <- recode(Subset$sex,"female" = "Female", "male"   = "Male")

plot_df <- Subset %>%
  filter(
    !is.na(ARG_div_shan), !is.na(sex), !is.na(age_category)) %>%
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
    !is.na(ARG_div_shan),
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

### 2.2 Linear model
With numeric age
```r
HIC_df <- plot_df %>% filter(Income_group == "HIC")
LMIC_df <- plot_df %>% filter(Income_group == "LMIC")

lm_HIC <- lm(ARG_div_shan ~ age_years + sex, data = HIC_df)
summary(lm_HIC)

lm_LMIC <- lm(ARG_div_shan ~ age_years + sex, data = LMIC_df)
summary(lm_LMIC)
```
Call:
lm(formula = ARG_div_shan ~ age_years + sex, data = HIC_df)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.73511 -0.28980  0.02956  0.32344  1.43123 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.6989980  0.0106261 159.889   <2e-16 ***
age_years    0.0058535  0.0002085  28.077   <2e-16 ***
sexMale     -0.0084761  0.0104798  -0.809    0.419    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4752 on 8254 degrees of freedom
Multiple R-squared:  0.08723,	Adjusted R-squared:  0.08701 
F-statistic: 394.4 on 2 and 8254 DF,  p-value: < 2.2e-16


Call:
lm(formula = ARG_div_shan ~ age_years + sex, data = LMIC_df)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.2185 -0.2862  0.0118  0.2661  1.3388 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.9532166  0.0279786  69.811  < 2e-16 ***
age_years    0.0037365  0.0005602   6.670 4.25e-11 ***
sexMale     -0.0276851  0.0265875  -1.041    0.298    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4174 on 984 degrees of freedom
Multiple R-squared:  0.04401,	Adjusted R-squared:  0.04207 
F-statistic: 22.65 on 2 and 984 DF,  p-value: 2.412e-10



### 2.3 Boxplot
With Categorised age
```r
plot_df <- Subset %>% filter(
    !is.na(ARG_div_shan),
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

ggplot(plot_df, aes(x = age_category, y = ARG_div_shan, fill = sex)) +
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
    title = "ARG Shannon Diversity by Age Category and Sex",
    x = "Age Category",
    y = expression("Shannon Diversity"),
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

### 2.4 LOESS
Change to GAM

--> Not statistically best, but nice for visualization.
```
ggplot(plot_df, aes(x = age_years, y = ARG_div_shan, color = sex)) +
  geom_point(alpha = 0.2, size = 1.0) +
  geom_smooth(method = "loess", se = TRUE, span = 0.6) +
  facet_wrap(~Income_group) +
  scale_color_npg() +
  labs(
    x = "Age (years)",
    y = "Shannon diversity",
    color = "Sex",
    title = "LOESS-smoothed ARG Shannon diversity by age, sex, and income group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA))
```

### 2.5 GAM 

```r
df_temp <- Subset |>
  dplyr::select(ARG_div_shan, age_years, sex, World_Bank_Income_Group) |>
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
    !is.na(Income_group))

plot_df <- Subset %>% filter(
    !is.na(ARG_div_shan),
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


fit <- gam(ARG_div_shan ~ sex + s(age_years, by = sex, k = 10) + sex*Income_group,
  method = "REML", data= df_temp)

summary(fit)
draw(fit)
```
Family: gaussian 
Link function: identity 

Formula:
ARG_div_shan ~ sex + s(age_years, by = sex, k = 10) + sex * Income_group

Parametric coefficients:
                           Estimate Std. Error t value
(Intercept)                1.913839   0.007415 258.112
sexfemale                  0.004505   0.010338   0.436
Income_groupLMIC           0.154708   0.022038   7.020
sexfemale:Income_groupLMIC 0.028337   0.031609   0.896
                           Pr(>|t|)    
(Intercept)                 < 2e-16 ***
sexfemale                     0.663    
Income_groupLMIC           2.37e-12 ***
sexfemale:Income_groupLMIC    0.370    
---
Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                         edf Ref.df      F p-value    
s(age_years):sexmale   4.795  5.839 110.76  <2e-16 ***
s(age_years):sexfemale 7.646  8.374  45.72  <2e-16 ***
---
Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.109   Deviance explained =   11%
-REML = 6076.7  Scale est. = 0.21631   n = 9244

### Draw the differences
```r
# 
diff <- difference_smooths(
  fit,
  smooth = "s(age_years)",
  n = 200)

draw(diff)
```

# Model check:
```r

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
  Income_group = levels(df_temp$Income_group))  # changed from World_Bank_Income_Group)

predictions <- predict(fit, newdata, se.fit = TRUE)
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit

diff_data <- newdata %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = c(pred, se))

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
    ARG_div_shan
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
    !is.na(ARG_div_shan),
    !is.na(Income_group),
    (sex == "Female" & age_years >= 15 & age_years <= 49) |
    (sex == "Male" & age_years >= 15 & age_years <= 49)
  )
```

### 4.2 Boxplot
```
ggplot(df_income_r_age_ARG, aes(x = Income_group, y = ARG_div_shan, fill = sex)) +
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

model_HIC <- lm(ARG_div_shan ~ age_years + sex,
  data = df_income_r_age_ARG %>% filter(Income_group == "HIC"))
summary(model_HIC)

model_LMIC <- lm(ARG_div_shan ~ age_years + sex,
  data = df_income_r_age_ARG %>% filter(Income_group == "LMIC"))
summary(model_LMIC)
```

## 2. ARG Load, numeric age, sex and income

### 2.1 Sample distribution
```r
Subset$sex <- recode(Subset$sex, "female" = "Female", "male"   = "Male")

plot_df <- Subset %>%
  filter(!is.na(ARG_div_shan), !is.na(sex), !is.na(age_years))

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
  filter(!is.na(ARG_div_shan),
         !is.na(sex),
         !is.na(age_years)) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))

ggplot(plot_df, aes(x = age_years, y = ARG_div_shan, color = sex, fill = sex)) +
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
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(age_years)) %>% nrow()

n_samples

lm_full <- lm(ARG_div_shan ~ sex + age_years, data = plot_df)
summary(lm_full)

```

### 2.6 Katan malli

# Differences 

```r
df_temp <- Subset %>%
  select(ARG_div_shan, sex, age_years, World_Bank_Income_Group) %>%
  filter(
    sex %in% c("male", "female"),
    !is.na(age_years), !is.na(World_Bank_Income_Group), !is.na(ARG_div_shan)) %>%
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
  ARG_div_shan ~ sex + s(age_years, by = sex, k = 10) + sex*Income_group,
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
plot(fit$fitted.values, df_temp$ARG_div_shan,
     xlab = "Fitted values",
     ylab = "Observed ARG_div_shan",
     main = "Response vs Fitted Values",
     pch = 16, col = "black")

# Add the diagonal line (y = x)
abline(a = 0, b = 1, col = "red", lwd = 2)

residuals <- df_temp$ARG_div_shan - fit$fitted.values
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
