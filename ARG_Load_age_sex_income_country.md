
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

## 2.Age 

### category sample distribution by continent 

```r
summary_table <- Subset %>%
  group_by(age_category, sex, country) %>%
  summarise(
    n_samples = n(),
    n_ARG_load = sum(!is.na(log10_ARG_load)),
    .groups = "drop"
  )
summary_table

summary_table <- Subset %>%
  count(continent, World_Bank_Income_Group, name = "n_samples")

```

summary_table <- Subset %>%
  group_by(continent, World_Bank_Income_Group) %>%
  summarise(
    n_samples = n(),
    sex_available = sum(!is.na(sex)),
    age_available = sum(!is.na(age_category)),
    ARG_available = sum(!is.na(log10_ARG_load)),
    .groups = "drop"
  )

  ```



