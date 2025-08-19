

###load necessary libraries###
pacman::p_load(tidyverse)



### Data Cleaning ###

### Step 1: Clean up birthplace variables to ensure consistency ###
#### Step 1a: Create new variable of yob and group the birthplace and residency into specific groups ###
# Clean birthplace into standard categories
survey_data <- read_csv("data/raw_data/AmericanSamoaDengueS_DATA_LABELS_2025-01-10_1548.csv") %>%
# survey_data <- read_csv("data/AmericanSamoaDengueS_DATA_LABELS_2024-02-27_1007.csv") %>%
      rename(birthplace=`4. Where was your child born?`,
             birthplace_other=`Where was your child born? Other:`) %>%
      mutate(birthplace_clean = case_when(
            birthplace == "American Samoa" ~ "American Samoa",
            is.na(birthplace_other) & birthplace == "American Samoa" ~ "American Samoa",  # Handle NA when birthplace is American Samoa
            str_detect(birthplace_other, "Samoa|samoa|Upolu|Saipipi") ~ "Western Samoa",
            str_detect(birthplace_other, "USA|Utah|Alaska|Washington|United States|California|Anchorage|San Francisco|Sacramento|WA|New York|GA") & 
                  !str_detect(birthplace_other, "Hawaii|HI|Honolulu|Oahu|Maui|Kauai|Big Island") ~ "USA",
            str_detect(birthplace_other, "Hawaii|HI|Honolulu|Oahu|Maui|Kauai|Big Island|Hawai'i|Hawaiian Islands|Hilo|Kona|Waipahu|Pearl Harbor") ~ "Hawaii",
            str_detect(birthplace_other, "Fiji") ~ "Fiji",
            str_detect(birthplace_other, "New Zealand") ~ "New Zealand",
            str_detect(birthplace_other, "Philippines|Phillipines") ~ "Philippines",
            str_detect(birthplace_other, "Australia") ~ "Australia",
            str_detect(birthplace_other, "Tuvalu") ~ "Tuvalu",
            TRUE ~ "Missing"
      )) %>%
          # mutate(birthplace_clean = case_when(
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("Western Samoa|Upolu|Saipipi|Tuanai Samoa|Savaii|Savai'i", ignore_case = TRUE)) ~ "Western Samoa",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("American Samoa", ignore_case = TRUE)) ~ "American Samoa",
          #   birthplace == "American Samoa" ~ "American Samoa",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("USA|United States|Utah|Alaska|Washington|California|Anchorage|San Francisco|Sacramento|New York|GA|NY", ignore_case = TRUE)) &
          #         !str_detect(birthplace_other, regex("Hawaii|HI|Honolulu|Oahu|Maui|Kauai|Big Island|Hilo|Kona|Waipahu|Pearl Harbor", ignore_case = TRUE)) ~ "USA",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("Hawaii|HI|Honolulu|Oahu|Maui|Kauai|Big Island|Hawai'i|Hawaiian Islands|Hilo|Kona|Waipahu|Pearl Harbor", ignore_case = TRUE)) ~ "Hawaii",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("Fiji", ignore_case = TRUE)) ~ "Fiji",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("New Zealand", ignore_case = TRUE)) ~ "New Zealand",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("Philippines|Phillipines", ignore_case = TRUE)) ~ "Philippines",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("Australia", ignore_case = TRUE)) ~ "Australia",
          #   !is.na(birthplace_other) & str_detect(birthplace_other, regex("Tuvalu", ignore_case = TRUE)) ~ "Tuvalu",
          #   TRUE ~ "Missing"
      #)) %>%
      mutate(yob = year(`Date of birth`)) %>%
      mutate(birth_group = case_when(
            birthplace_clean == "Missing" ~ NA_character_,
            birthplace_clean == "American Samoa" ~ "American Samoa",
            # birthplace_clean %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            # birthplace_clean %in% c("New Zealand", "USA", "Australia") ~ "New Zealand, USA and Australia",
            birthplace_clean %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines") ~ "Pacific Islands",
            birthplace_clean %in% c("New Zealand", "USA", "Australia","Hawaii") ~ "New Zealand, USA and Australia",
            TRUE ~ "Other"
      ))#%>%
      
      # Extract year of birth
      # mutate(yob = year(`Date of birth`)) %>%
      
      # Assign to birth groups
      # mutate(birth_group = case_when(
      #       birthplace_clean == "Missing" ~ NA_character_,
      #       birthplace_clean == "American Samoa" ~ "American Samoa",
      #       birthplace_clean %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
      #       birthplace_clean %in% c("New Zealand", "USA", "Australia") ~ "New Zealand, USA and Australia",
      #       TRUE ~ "Other"
      # ))

###Step 1b:Create a new dataframe with cleaned up residency variables###
survey_data1 <- survey_data %>%
###Step 1c: Remove the indeterminate cases###
  filter (`Final Test Result`!= "I") %>%
### Step 1d:Remove the following criteria 1) Missing birthplace 2) Missing residency###
filter(`Record ID` != "477") %>%
filter(`Record ID` != "763") %>%
filter(`Record ID` != "135") %>%
### Step 1e: Clean up the residency variables###
  mutate(
    `5. My child has been a resident of American Samoa` = case_when(
      `Record ID` == 739 & `5. My child has been a resident of American Samoa` == "Form is blank" ~ "All his/her life", 
      TRUE ~`5. My child has been a resident of American Samoa`
    ),
    #category is ignoring residency history and including location information
    #replacing initial AS residency marked as birth year, and include location data
    `My child has been a resident of American Samoa since (year)` = case_when(
      `Record ID` == 279 & `My child has been a resident of American Samoa since (year)` == 2011 ~ 2023, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `5. My child has been a resident of American Samoa` = case_when (
      `Record ID` == 279 & `5. My child has been a resident of American Samoa` == "Since __ 2011 __ (year)" ~"Since __ 2023 __ (year)", 
      TRUE ~`5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when (
      `Record ID` == 342 & `My child has been a resident of American Samoa since (year)` == 2010~2017, 
      TRUE ~`My child has been a resident of American Samoa since (year)`
    ),
    `5. My child has been a resident of American Samoa` = case_when (
      `Record ID` == 342 & `5. My child has been a resident of American Samoa` == "Since __ 2010 __ (year)" ~"Since __ 2017 __ (year)", 
      TRUE ~`5. My child has been a resident of American Samoa`
    ),
      `5. My child has been a resident of American Samoa` = case_when(
        `Record ID` == 415 & `5. My child has been a resident of American Samoa` == "All his/her life" ~ "", 
        TRUE ~ `5. My child has been a resident of American Samoa`
    ),
    `5. My child has been a resident of American Samoa` = case_when(
      `Record ID` == 322 & `5. My child has been a resident of American Samoa` == "All his/her life" ~ "Since __ 2020 __ (year)", 
      TRUE ~ `5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when(
      `Record ID` == 322 & is.na(`My child has been a resident of American Samoa since (year)`) ~ 2020, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `5. My child has been a resident of American Samoa` = case_when(
      `Record ID` == 487 & `5. My child has been a resident of American Samoa` == "All his/her life" ~ "Since __ 2023 __ (year)", 
      TRUE ~ `5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when (
      `Record ID` == 487 & is.na(`My child has been a resident of American Samoa since (year)`) ~2023, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 487 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 487 & is.na(`Residence in [residence_other2], year start`)~ 2010, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 487 & is.na(`Residence in [residence_other2], year end`) ~ 2018, 
      TRUE ~`Residence in [residence_other2], year end`
    ),
    `Residence in [residence_other], year start`= case_when (
      `Record ID` == 247 & `Residence in [residence_other], year start` == 2011 ~ 2009, 
      TRUE ~ `Residence in [residence_other], year start`
    ),
    `Residence in [residence_other], year start` = case_when(
      `Record ID` == 799 & `Residence in [residence_other], year start` == 2013 ~ 2007, 
      TRUE ~ `Residence in [residence_other], year start`
    ),
    `5. My child has been a resident of American Samoa` = case_when(
      `Record ID` == 434 & `5. My child has been a resident of American Samoa` == "Since __ 2009 __ (year)" ~ "Since __ 2023 __ (year)", 
      TRUE ~ `5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when (
      `Record ID` == 434 & `My child has been a resident of American Samoa since (year)` == 2009 ~ 2023, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 434 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 434 & is.na(`Residence in [residence_other2], year start`) ~ 2009, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 434 & is.na(`Residence in [residence_other2], year end`) ~ 2018, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 348 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 348 & is.na(`Residence in [residence_other2], year start`) ~ 2010, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 348 & is.na(`Residence in [residence_other2], year end`) ~ 2011, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `If your child has lived somewhere else before becoming a resident of American Samoa, where did they live before?` = case_when (
      `Record ID` == 371 & is.na(`If your child has lived somewhere else before becoming a resident of American Samoa, where did they live before?`) ~ "Western Samoa", 
      TRUE ~`If your child has lived somewhere else before becoming a resident of American Samoa, where did they live before?`
    ),
    `Residence in [residence_other], year start` = case_when (
      `Record ID` == 371 & is.na(`Residence in [residence_other], year start`) ~ 2011, 
      TRUE ~ `Residence in [residence_other], year start`
    ),
    `Residence in [residence_other], year end` = case_when (
      `Record ID` == 371 & is.na(`Residence in [residence_other], year end`) ~ 2019, 
      TRUE ~ `Residence in [residence_other], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 480 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 480 & is.na(`Residence in [residence_other2], year start`) ~ 2011, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 480 & is.na(`Residence in [residence_other2], year end`) ~ 2016, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `5. My child has been a resident of American Samoa` = case_when(
      `Record ID` == 213 & `5. My child has been a resident of American Samoa` == "Since __ 2014 __ (year)" ~ "Since __ 2023 __ (year)", 
      TRUE ~ `5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when (
      `Record ID` == 213 & `My child has been a resident of American Samoa since (year)` == 2014 ~ 2023, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 213 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 213 & is.na(`Residence in [residence_other2], year start`)~ 2014, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 213 & is.na(`Residence in [residence_other2], year end`) ~ 2018, 
      TRUE ~`Residence in [residence_other2], year end`
    ),    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...44` = case_when (
      `Record ID` == 542 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...44`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...44`
    ),
    `Residence in [residence_other3], year start` = case_when (
      `Record ID` == 542 & is.na(`Residence in [residence_other3], year start`) ~ 2011, 
      TRUE ~ `Residence in [residence_other3], year start`
    ),
    `Residence in [residence_other3], year end` = case_when (
      `Record ID` == 542 & is.na(`Residence in [residence_other3], year end`) ~ 2012, 
      TRUE ~ `Residence in [residence_other3], year end`
    ),
    `Residence in [residence_other], year end` = case_when (
      `Record ID` == 542 & `Residence in [residence_other], year end` == 2014~2015,
      TRUE ~ `Residence in [residence_other], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 608 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 608 & is.na(`Residence in [residence_other2], year start`) ~ 2013, 
      TRUE~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when(
      `Record ID` == 608 & is.na(`Residence in [residence_other2], year end`) ~ 2015, 
      TRUE~ `Residence in [residence_other2], year end`
      ),
    `Residence in [residence_other3], year start` = case_when (
      `Record ID` == 640 & `Residence in [residence_other3], year start` == 2019 ~ 2009, 
      TRUE ~`Residence in [residence_other3], year start`
    ),
    `Residence in [residence_other3], year end` = case_when (
      `Record ID` == 640 & is.na(`Residence in [residence_other3], year end`) ~ 2010, 
      TRUE ~ `Residence in [residence_other3], year end`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 685 & is.na(`Residence in [residence_other2], year end`) ~ 2023, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `5. My child has been a resident of American Samoa` = case_when (
      `Record ID` == 686 & `5. My child has been a resident of American Samoa` == "All his/her life"~ "Since __ 2023 __ (year)", 
      TRUE ~ `5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when (
      `Record ID` == 686 & is.na(`My child has been a resident of American Samoa since (year)`) ~ 2023, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 686 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 686 & is.na(`Residence in [residence_other2], year start`) ~ 2014, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 686 & is.na(`Residence in [residence_other2], year end`) ~ 2021, 
      TRUE ~`Residence in [residence_other2], year end` 
    ),
    `5. My child has been a resident of American Samoa` = case_when(
      `Record ID` == 696 & `5. My child has been a resident of American Samoa` == "All his/her life" ~ "Since __ 2023 __ (year)", 
      TRUE ~  `5. My child has been a resident of American Samoa`
    ),
    `My child has been a resident of American Samoa since (year)` = case_when (
      `Record ID` == 696 & is.na(`My child has been a resident of American Samoa since (year)`) ~ 2023, 
      TRUE ~ `My child has been a resident of American Samoa since (year)`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 696 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 696 & is.na(`Residence in [residence_other2], year start`) ~ 2010, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 696 & is.na(`Residence in [residence_other2], year end`) ~ 2017,
      TRUE ~   `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 718 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when(
      `Record ID` == 718 & is.na(`Residence in [residence_other2], year start`) ~ 2009, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 718 & is.na(`Residence in [residence_other2], year end`) ~ 2012, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `Residence in [residence_other], year start` = case_when (
      `Record ID` == 722 & `Residence in [residence_other], year start` == 2015 ~ 2008, 
      TRUE ~ `Residence in [residence_other], year start`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 738 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 738 & is.na(`Residence in [residence_other2], year start`) ~  2009, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 738 & is.na(`Residence in [residence_other2], year end`) ~ 2017, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 741 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when(
      `Record ID` == 741 & is.na(`Residence in [residence_other2], year start`) ~ 2009, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when(
      `Record ID` == 741 & is.na(`Residence in [residence_other2], year end`) ~ 2012, 
      TRUE ~`Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 747 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 747 & is.na(`Residence in [residence_other2], year start`) ~ 2008, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 747 & is.na(`Residence in [residence_other2], year end`) ~ 2010, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 771 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 771 & is.na(`Residence in [residence_other2], year start`) ~ 2008, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 771 & is.na(`Residence in [residence_other2], year end`) ~ 2015, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...44` = ifelse (
      `Record ID` == 853 & `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...44` == "Colorado",
      "American Samoa", 
      `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...44`
    ),
    `Residence in [residence_other3], year start` = ifelse(
      `Record ID` == 853 & `Residence in [residence_other3], year start`== 2014, 
      2008, 
      `Residence in [residence_other3], year start`
    ),
    `Residence in [residence_other3], year end` = ifelse(
      `Record ID` == 853 & `Residence in [residence_other3], year end` == 2019, 
      2013, 
      `Residence in [residence_other3], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = ifelse (
      `Record ID` == 853 & `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` == "Hawaii",
      "Hawaii", 
      `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = ifelse(
      `Record ID` == 853 & `Residence in [residence_other2], year start`== 2019, 
      2014, 
      `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = ifelse(
      `Record ID` == 853 & `Residence in [residence_other2], year end`== 2020, 
      2020, 
      `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 800 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~"American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 800 & is.na(`Residence in [residence_other2], year start`) ~ 2019, 
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 800 & is.na(`Residence in [residence_other2], year end`) ~ 2023, 
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41` = case_when (
      `Record ID` == 719 & is.na(`If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`) ~ "American Samoa", 
      TRUE ~ `If your child has live somewhere else before becoming a resident of American Samoa, where did they live before?...41`
    ),
    `Residence in [residence_other2], year start` = case_when (
      `Record ID` == 719 & is.na(`Residence in [residence_other2], year start`) ~ 2009,
      TRUE ~ `Residence in [residence_other2], year start`
    ),
    `Residence in [residence_other2], year end` = case_when (
      `Record ID` == 719 & is.na(`Residence in [residence_other2], year end`) ~ 2023,
      TRUE ~ `Residence in [residence_other2], year end`
    ),
    `5. My child has been a resident of American Samoa` = case_when (
      `Record ID` == 719 & `5. My child has been a resident of American Samoa` == "All his/her life" ~ "", 
      TRUE ~`5. My child has been a resident of American Samoa`
    ),
    `5. My child has been a resident of American Samoa` = case_when (
      `Record ID` == 895 & `5. My child has been a resident of American Samoa` == "All his/her life" ~ "", 
      TRUE ~ `5. My child has been a resident of American Samoa`
    )#,
    #   birthplace_clean = case_when(
    #     `Record ID` == 135 & birthplace_clean == "American Samoa" ~ "Western Samoa",
    #     TRUE ~ birthplace_clean
    #   ),
    #   `Residence in [residence_other], year start` = case_when(
    #     `Record ID` == 135 & `Residence in [residence_other], year start` == 2007 ~ 2010,
    #     TRUE ~ `Residence in [residence_other], year start`
    # )
    )



###Step 2: Create location matrices with start and end dates based on residency info ###
survey_data1 <- survey_data1 %>%
  mutate(across(c(38:40, 41:43, 44:46)))

# Create location_matrix_1
location_matrix_1 <- survey_data1 [, c(1, 3, 4, 38:40)] %>%
filter(`ID Color` != "Mop-up (Not tested yet)")
colnames(location_matrix_1) <- c("Record ID", "ID Color", "ID Number", "location", "start", "end")

# Create location_matrix_2
location_matrix_2 <- survey_data1[, c(1, 3, 4, 41:43)] %>%
  filter(`ID Color` != "Mop-up (Not tested yet)")
colnames(location_matrix_2) <- c("Record ID", "ID Color", "ID Number", "location", "start", "end")

# Create location_matrix_3
location_matrix_3 <- survey_data1[, c(1, 3, 4, 44:46)] %>%
  filter(`ID Color` != "Mop-up (Not tested yet)")
colnames(location_matrix_3) <- c("Record ID", "ID Color", "ID Number", "location", "start", "end")

# Combine the matrices
location_matrix <- bind_rows(location_matrix_1, location_matrix_2, location_matrix_3) %>%
  filter(!is.na(location))


### Step 3: Clean up the locations for the combined location matrix ###
clean_location <- location_matrix %>%
  rename(location_raw = location) %>%
  mutate(location = case_when(
    str_detect(location_raw, "American Samoa") ~ "American Samoa",
    str_detect(location_raw, "Samoa|samoa|Savaii") & !str_detect(location_raw, "American Samoa") ~ "Western Samoa",
    str_detect(location_raw, "USA|Alaska|Washington|United States|California|Colorado|UTAH") & 
      !str_detect(location_raw, "Hawaii|HI|Honolulu") ~ "USA",  # Exclude Hawaii mentions from USA
    str_detect(location_raw, "Honolulu|Hawaii|HI|Oahu|Maui|Kauai|Big Island|Hawaiian Islands|Hawai'i") ~ "Hawaii",
    str_detect(location_raw, "New Zealand") ~ "New Zealand",
    str_detect(location_raw, "Philippines") ~ "Philippines",
    str_detect(location_raw, "Australia") ~ "Australia",
    str_detect(location_raw, "Tuvalu") ~ "Tuvalu",
    str_detect(location_raw, "Fiji") ~ "Fiji",
    TRUE ~ "Missing"
  ))

### Step 4: Combining all residency information ###
base <- expand.grid(`Record ID`=survey_data$`Record ID`,year=2000:2023)

###Step 4a: Create a dataframe for all possibilities with year and id combinations###
possibilities <- base %>% left_join(survey_data1, by = "Record ID")%>%
  filter(year>=yob) %>%
  select(`Record ID`,yob, birthplace_clean, year) 

###Step 4b: Create a dataframe with place of birth and year of birth from cleaned dataset#
pob1 <-select(survey_data1, `Record ID`, year=yob, location=birthplace_clean)


###Step 4c: Create a dataframe which shows individuals who have lived all their lives in American Samoa#
asresidency1 <- base %>% left_join (survey_data1, by = "Record ID")%>%
  filter(year>yob & `5. My child has been a resident of American Samoa` == "All his/her life" )%>%
  select(`Record ID`,year) %>%
  mutate(location = "American Samoa")

### Step 4d: Merge the data frame for who have resided in American Samoa from a specific year#
asresidency2 <- base %>% left_join (survey_data1, by = "Record ID")%>%
  filter(!is.na(`My child has been a resident of American Samoa since (year)`))%>%
  filter(year>=`My child has been a resident of American Samoa since (year)`) %>%
  select(`Record ID`,year) %>%
  mutate(location = "American Samoa")


### Step 4e: Join the datasets with clean_locations and base###
### Ignore warning message###
lom3 <- base %>% left_join(clean_location, by = "Record ID")%>%
filter(year>=start,year<=end) %>%
select(`Record ID`, year, location) 

### Step 4f: Each individual assigned to be living in AS in 2023 ###
 as2023 <-select(survey_data1, `Record ID`)%>%
  mutate(year=2023, location= "American Samoa")

 
### Step 4g: combining the datsets for place of birth, location matrix, as residency all their life and since specific year# 
 combination <- bind_rows(pob1, lom3, asresidency1, asresidency2, as2023)%>%
   distinct()
 
 
 combination2 <- possibilities %>% left_join(combination, by = c("Record ID", "year"))%>%
   arrange(`Record ID`, year)%>%
   group_by(`Record ID`, year)%>%
   mutate(n=n()) %>%
   mutate(location= paste(location,collapse=", "))%>%
   mutate(location = ifelse(location=="NA", NA, location))%>%
   distinct()
 
 ###Step 5: Separate the location into distinct locations either former or latter location as some individuals lived in 1 or more places in one year ###
 combination2_sep <- combination2 %>%
   separate(location, into = c("Former_location", "Latter_location"),sep = ", ", fill = "right") %>%
   mutate(Latter_location = ifelse(is.na(Latter_location) & !is.na(Former_location), Former_location, Latter_location)) %>%
   distinct() %>%
   ungroup()
 
 
 #create long dataset with former and latter location and groups 
long_dataset <-  combination2_sep %>%
                        select(`Record ID`, yob, birthplace_clean, year, Former_location,Latter_location)%>%
      mutate(Former_location_group = case_when(
            Former_location == "Missing" ~ NA_character_,
            Former_location == "American Samoa" ~ "American Samoa",
            # Former_location %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            # Former_location %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
            Former_location %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines") ~ "Pacific Islands",
            Former_location %in% c("New Zealand", "USA", "Australia","Hawaii") ~ "Non-endemic",
            TRUE ~ Former_location
      ))%>%
      mutate(Latter_location_group = case_when(
            Latter_location == "Missing" ~ NA_character_,
            Latter_location == "American Samoa" ~ "American Samoa",
            # Latter_location %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            # Latter_location %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
            Latter_location %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines") ~ "Pacific Islands",
            Latter_location %in% c("New Zealand", "USA", "Australia", "Hawaii") ~ "Non-endemic",
            TRUE ~ Latter_location
      ))

write.csv(long_dataset, "data/generated_data/long_dataset.csv")

#select the variables interested in from survey data
wide_dataset <- survey_data1 %>%
      select("Record ID", "Date of birth", "Calculated age based on DOB.", 
             "School", "Sex", "Village of Residence", "Final Test Result", "birthplace_clean","birth_group","yob",
             "5. My child has been a resident of American Samoa"
             )

write.csv(wide_dataset, "data/generated_data/wide_dataset.csv")
      
 
 
 
# create first dataset with former location
#  former_location_data <- combination2_sep %>%
#    select(`Record ID`, yob, birthplace_clean, year, Former_location)
#  
#  #create second dataset with latter location
#  latter_location_data <- combination2_sep %>%
#    select(`Record ID`, yob, birthplace_clean, year, Latter_location)
#  
#  # Add location grouping to latter_location_data
#  latter_long_data <- latter_location_data %>%
#    mutate(location_group = case_when(
#      Latter_location == "Missing" ~ NA_character_,
#      Latter_location == "American Samoa" ~ "American Samoa",
#      Latter_location %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
#      Latter_location %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
#      TRUE ~ Latter_location
#    ))
#  
#  # Create wide format for locations
#  latter_location_wider <- latter_long_data %>%
#    select(`Record ID`, yob, birthplace_clean, year, Latter_location) %>%
#    pivot_wider(names_from = year, values_from = Latter_location)
#  
#  # Create wide format for location groups
#  latter_location_groups_wider <- latter_long_data %>%
#    select(`Record ID`, year, location_group) %>%
#    pivot_wider(
#      names_from = year,
#      values_from = location_group,
#      names_prefix = "location_group_"
#    )
#  
#  # Pivot wider for former location
#  former_location_wider <- former_location_data %>%
#    pivot_wider(names_from = year, values_from = Former_location)
#  
#  #select the variables interested in from survey data
#  survey_data2 <- survey_data1 %>%
#    select("Record ID", "Date of birth", "Calculated age based on DOB.", "School", "Sex", "Village of Residence", "Final Test Result", "birth_group")
#  
#  # Create merged datasets
#  merged_former1 <- inner_join(former_location_wider, survey_data2, by = "Record ID") 
#  merged_latter1 <- inner_join(latter_location_wider, survey_data2, by = "Record ID")
#  
#  # Create final analytical dataset with both locations and location groups
#  analytical_dataset <- survey_data2 %>%
#    left_join(
#      merged_latter1 %>% 
#        select(`Record ID`, 
#               yob, 
#               birthplace_clean,
#               starts_with("20")),  # Original location columns
#      by = "Record ID"
#    ) %>%
#    left_join(
#      latter_location_groups_wider,  # Add location group columns
#      by = "Record ID"
#    )
#  
# 
# 
# #export datasets
# write.csv(merged_former1, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\merged_former1.csv")
# #export datasets
# write.csv(merged_latter1, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\merged_latter1.csv")
# 
# write.csv(combination2, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\combination2.csv")
# 
# write.csv(survey_data2, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\survey_data1.csv")
# 
# write.csv(latter_location_data, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\latter_location_data.csv")
# 
# write.csv(latter_long_data, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\latter_long_data.csv")
# 
#write.csv(analytical_dataset, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\analytical_dataset.csv")

