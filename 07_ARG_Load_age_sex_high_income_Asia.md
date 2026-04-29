
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

```
--> 330 samples removed

## 2.High income countries in EU

### 2.1 Sample distribution

--> Filter China it is labeles as Upper-Middle income country
```r
hic_As <- Subset %>%
  filter(
    World_Bank_Income_Group == "High income",
    continent == "Asia")

unique(hic_As$country)

```
Age category (for boxplot)

```r
hic_As_clean <- hic_As %>%
  filter(
    !is.na(age_category),
    !is.na(sex))
table(hic_As_clean$sex, useNA = "ifany")

table_df <- hic_As_clean %>%
  count(age_category, sex) %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(Total = female + male)

table_df <- table_df %>%
  mutate(
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"
    ))
  ) %>% arrange(age_category)

table_image <- table_df %>%
  gt() %>%
  cols_label(
    age_category = "Age Category",
    female = "Female (n)",
    male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(
    title = "High-Income Asia: Sample Distribution by Age and Sex"
  )
```
Numeric age (for linear regression)

```r
hic_As_clean_num <- hic_As %>% filter(!is.na(age_years), !is.na(sex))

sex_table <- hic_As_clean_num %>% count(sex)

sex_gt <- sex_table %>%
  gt() %>%
  cols_label(
    sex = "Sex",
    n = "Count (n)"
  ) %>%
  tab_header(title = "Sex distribution in the High-Income Asian cohort")
```


## 2.2 Boxplot

```r
plot_df <- hic_As_clean %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(age_category)
  ) %>%
  mutate(
    sex = case_when(
      sex %in% c("female", "Female") ~ "Female",
      sex %in% c("male", "Male") ~ "Male",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"
    ))
  ) %>%
  filter(!is.na(sex))

ggplot(plot_df, aes(x = age_category, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by Age Category and Sex (High-Income Asia)",
    x = "Age Category",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  stat_compare_means(
    aes(group = sex),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

```

### 2.3 Linear regression
```r
model <- lm(
  log10_ARG_load ~ sex + age_years,
  data = hic_As_clean_num)

summary(model)
```

