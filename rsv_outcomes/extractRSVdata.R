library("tidyr")
library("dplyr")

# load data from RSV Model output
file_path <- "./rsv_outcomes/outcomes_all.RDS"
load(file_path)

age_0to5mo <- c("<1mo", "1mo", "2mo", "3mo", "4mo", "5mo")
age_6to23mo <- c("6mo", "7mo", "8mo", "9mo", "10mo", "11mo", "1yr")
age_2to4yr <- c("2yr", "3yr", "4yr")
age_5to17y <- c("5-9yr", "10-14yr", "15-17yr")
age_18to100y <- c("18-24yr", "25-34yr", "35-44yr", "45-54yr", "55-64yr", "65-74yr", "75+yr")


# grab the annual GP visits data
GPvisits <- df_outcomes$health_outcomes %>% 
  filter(effects == "total", metric == "GP visits") %>%
  select(seed, age_group, value, inter) %>%
  mutate(inter = paste(inter))

GPvisits_table <- GPvisits %>%
  pivot_wider(names_from = age_group, values_from = value) %>%
  mutate(`15-17yr` = `15-24yr` * (3/10)) %>% # add column for 15-17yo
  mutate(`18-24yr` = `15-24yr` * (7/10)) %>% # add column for 18-24yo
  select(c(seed:`10-14yr`,`15-17yr`, `18-24yr`, `25-34yr`, `35-44yr`, `45-54yr`, `55-64yr`, `65-74yr`, `75+yr`)) %>% # grab only those columns needed
  mutate(age_0to5mo = rowSums(across(age_0to5mo))) %>%
  mutate(age_6to23mo = rowSums(across(age_6to23mo))) %>%
  mutate(age_2to4yr = rowSums(across(age_2to4yr))) %>%
  mutate(age_5to17y = rowSums(across(age_5to17y))) %>%
  mutate(age_18to100y = rowSums(across(age_18to100y))) %>%
  select(seed, inter, age_0to5mo, age_6to23mo, age_2to4yr, age_5to17y, age_18to100y) %>%
  pivot_longer(col = !(seed | inter), names_to = "age_group", values_to = "cases") %>%
  pivot_wider(names_from = inter, values_from = cases) # make a table with interventions as columns
  

# summary table
# GPvisits_summary <- GPvisits_table %>%
#   group_by(inter, age_group) %>%
#   summarise(GPvisits_mean = mean(cases),
#             GPvisits_loCI = my_quantile(cases, 0.025),
#             GPvisits_hiCI = my_quantile(cases, 0.975))


# calculate GP visits averted
intervention_list <- paste(unique(df_outcomes$health_outcomes[["inter"]]))
new_interventions <- intervention_list[!intervention_list %in% c("PAL_VHR_S","VAC_65_S","VAC_75_S")]
status_quo <- intervention_list[1]
status_quo_col <- unlist(unname(GPvisits_table[status_quo]))


# fn to calculate the averted number of cases
avert_fun <- function(intervention, status_quo = status_quo){
  (status_quo - intervention) / status_quo
}


# calculate the GP visits averted using the avert fn above (as samples)
GPvisits_averted <- GPvisits_table %>% 
  mutate(across(new_interventions, ~avert_fun(.x, status_quo = status_quo_col))) %>%
  # mutate(status_quo_GP_visits = PAL_VHR_S)  %>%
  select(-c(PAL_VHR_S, VAC_65_S, VAC_75_S)) %>% 
  pivot_longer(col = !(seed | age_group), names_to = "intervention", values_to = "prop_averted") 

# summary table
GPvisits_averted_summary <- GPvisits_averted %>%
  group_by(intervention, age_group) %>%
  summarise(GPvisits_prop_averted_mean = mean(prop_averted),
            GPvisits_prop_averted_loCI = my_quantile(prop_averted, 0.025),
            GPvisits_propaverted_hiCI = my_quantile(prop_averted, 0.975))
  

# saveRDS(GPvisits_averted, "data_ProportionGPVisitsAverted.RDS")
