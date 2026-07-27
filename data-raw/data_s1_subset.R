## code to prepare `data_s1_subset` dataset goes here
library(dplyr)

# 1. Read the data
# Ideally, move 'Data.csv' into the 'data-raw' folder so this path is relative
# If the file is in 'data-raw/', use:
# raw_path <- "data-raw/Data.csv"
# For now, using your specific path:
raw_path <- '/Users/cliattx/Downloads/Data.csv'

if(!file.exists(raw_path)) stop("Raw data file not found!")

data_demo <- read.csv(raw_path)

# 2. Clean the data (Your exact code)
data_demo <- data_demo %>%
  filter(Outcome == 'DCPM') %>%
  mutate(Session = as.numeric(Session)) %>%
  group_by(Subject) %>%
  mutate(
    First_Treatment_Session = min(Session[Phase == 1], na.rm = TRUE),
    Time_Centered = Session - First_Treatment_Session + 1,
    Time.c_phase = Time_Centered * Phase
  ) %>%
  mutate(
    Time.c_phase = ifelse(Time.c_phase == 0, 0, Time.c_phase - 1)
  ) %>%
  ungroup()

data_s1_subset <- data_demo %>%
  select(Subject, Session, Y = Score, Phase, Time = Time_Centered, Time.c_phase)

data_s1_subset <- as.data.frame(data_s1_subset)

# 3. Save the data to the package
usethis::use_data(data_s1_subset, overwrite = TRUE)
