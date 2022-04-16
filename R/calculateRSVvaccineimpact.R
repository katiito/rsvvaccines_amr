##########################################################
# calculate the impact of RSV vaccination on prescriptions and antibiotic prescribing in the UK
# Code for the Paper: Atkins, Hodgson, Jit and Davies (2021) medRxiv
##########################################################

###########################################################
### 0. Set up libraries and parameters for use later
###########################################################

# load some libraries
library(here)
library(tidyverse)
library(gridExtra)
library(data.table)

# read in RSV model output and required functions
source(here::here("R/RSV_AMR_functions.R")) # TAKEN / ADAPTED FROM CHAE et al. (2018) Emerging Infectious Diseases
source(here::here("rsv_outcomes/extractRSVdata.R")) # TAKEN FROM HODGSON et al. (2020) BMC Med

# number of samples used to generate the uncertainty intervals
n <- 1000

# GET POPULATION SIZES - England and Scotland from 2019 ONS website (saved as data/ONS_England_agedist.csv)
# ENGLAND
engpopsize0to5mo <- 618858 * 0.5 # assumption that 0-1y is equally split between <6mo and >6mo
engpopsize6to23mo <- 618858 * 0.5 + 644056 # assumption that 0-1y is equally split between <6mo and >6mo
engpopsize2to4yr <- 2036723
engpopsize5to17yr <- 8723931
engpopsize18to100yr <- 44263393
engpopsize_0to4y <- engpopsize0to5mo + engpopsize6to23mo + engpopsize2to4yr
engpopsize0to100yr <- engpopsize0to5mo + engpopsize6to23mo + engpopsize2to4yr + engpopsize5to17yr + engpopsize18to100yr
ew_popsize0to100yr_2018 <- 66436000

popsize <- list(engpopsize0to5mo = engpopsize0to5mo,
                engpopsize6to23mo = engpopsize6to23mo,
                engpopsize2to4yr = engpopsize2to4yr,
                engpopsize5to17yr = engpopsize5to17yr,
                engpopsize18to100yr = engpopsize18to100yr,
                engpopsize0to100yr = engpopsize0to100yr)

# SCOTLAND
scotpopsize0to5mo <- 50772 * 0.5
scotpopsize6to11mo <- 50772 * 0.5
scotpopsize12to23mo <- 52734
scotpopsize2to4yr <- 168209
scotpopsize5to17yr <- 757447
scotpopsize18to100yr <- 4434138

############################################################### 
#### 1. Calculate the proportion of GP visits averted for  #### 
####          each vaccine strategy for each age group     #### 
###############################################################

# sample of annual proportion of GP visits averted by age group, columns by vaccine strategy

age_groups <- GPvisits_averted$age_group %>% unique # Get the age groups names
# Generate the tables
gppropaverted <- age_groups %>% map(
  ~(GPvisits_averted %>%
    filter(age_group == .x) %>%
    pivot_wider(names_from = intervention, values_from = prop_averted) %>%  
    select(!(age_group | seed))
  )
)

names(gppropaverted) <- c(
  "GPvisits_proportionaverted_0to5mo",
  "GPvisits_proportionaverted_6to23mo",
  "GPvisits_proportionaverted_2to4yr",
  "GPvisits_proportionaverted_5to17yr",
  "GPvisits_proportionaverted_18to100yr")

# PRINT TABLE 1 FROM PAPER
PrintAgeStratified(GPvisits_averted_summary)


##################################################################################  
#### 2. i. Calculate the number of GP visits that result in an                #### 
###     ABX prescription (any class) attributable to RSV for each age group   #### 
##################################################################################

# draw n samples from the annual RSV-attributable Rxs (per 100,000 persons) for each age group 
# under BASE CASE assumptions i.e. from Taylor et al. assuming CIs indicative of 95% CIs of a triangular distribution

Rx_0to5mo_per100000py_a <- rtriang(n = n, mode = 8328, l = 5547, r = 10265)
Rx_6to23mo_per100000py_a <- rtriang(n = n, mode = 11916 , l = 8432, r = 13684)
Rx_2to4y_per100000py_a <- rtriang(n = n, mode = 7495, l = 5084, r = 9051)
Rx_5to17y_per100000py_a <- rtriang(n = n, mode = 1091, l = 686, r = 1427)

# Taylor et al. don't calculate the RSV-attributable Rxs for 18+y so we assume:
# i) the rate of Rx is is *same* as that for 5-17y
Rx_18to100y_per100000py_a <- rtriang(n = n, mode = 1091, l = 686, r = 1427)

rx_rates_A <- list(Rx_0to5mo_per100000py = Rx_0to5mo_per100000py_a,
                   Rx_6to23mo_per100000py = Rx_6to23mo_per100000py_a,
                   Rx_2to4y_per100000py = Rx_2to4y_per100000py_a,
                   Rx_5to17y_per100000py = Rx_5to17y_per100000py_a,
                   Rx_18to100y_per100000py = Rx_18to100y_per100000py_a)


###########################################################
## 2. ii. Calculate the Number of Rx associated with Abx (*_b denotes Sensivitiy Analysis I, _c denotes Sensivity Analysis II)
### Fitzpatrick 2020 (CID) Scottish Child study https://doi.org/10.1093/cid/ciaa403
###########################################################
Rx_0to11mo_all_per100000y <- 532.7 * 100
Rx_1to4y_all_per100000y <- 628.5 * 100

pc_RxRSV_0to11mo <- rtriang(n = n, mode = 0.052, l = 0.039, r = 0.064)
pc_RxRSV_1to4y <- rtriang(n = n, mode = 0.058, l = 0.046, r = 0.070) 

Rx_0to11mo_per100000y_b <- Rx_0to11mo_all_per100000y * pc_RxRSV_0to11mo
Rx_1to4y_per100000y_b <- Rx_1to4y_all_per100000y * pc_RxRSV_1to4y

# Now calculate the finer age stratification rates
# - assume rates within larger age groups are the same.
# - assume outside-study rates are proportionally the same as the English study (ie ratio of 5-17y to 2-4y Rx rates are the same across studies)

Rx_0to5mo_per100000py_b <- Rx_0to11mo_per100000y_b
Rx_6to23mo_per100000py_b <- Rx_0to11mo_per100000y_b * (scotpopsize6to11mo / (scotpopsize6to11mo + scotpopsize12to23mo)) +
                          Rx_1to4y_per100000y_b * (scotpopsize12to23mo / (scotpopsize6to11mo + scotpopsize12to23mo)) 
Rx_2to4y_per100000py_b <- Rx_1to4y_per100000y_b
Rx_5to17y_per100000py_b <- (Rx_5to17y_per100000py_a / Rx_2to4y_per100000py_a) * Rx_2to4y_per100000py_b
Rx_18to100y_per100000py_b <- Rx_5to17y_per100000py_b

rx_rates_B <- list(Rx_0to5mo_per100000py = Rx_0to5mo_per100000py_b,
                   Rx_6to23mo_per100000py = Rx_6to23mo_per100000py_b,
                   Rx_2to4y_per100000py = Rx_2to4y_per100000py_b,
                   Rx_5to17y_per100000py = Rx_5to17y_per100000py_b,
                   Rx_18to100y_per100000py = Rx_18to100y_per100000py_b)

rx_rates_C <- list(Rx_0to5mo_per100000py = Rx_0to5mo_per100000py_b,
                   Rx_6to23mo_per100000py = Rx_6to23mo_per100000py_b,
                   Rx_2to4y_per100000py = Rx_2to4y_per100000py_b,
                   Rx_5to17y_per100000py = Rx_5to17y_per100000py_b,
                   Rx_18to100y_per100000py = 0 * Rx_18to100y_per100000py_b) #re-scale so adults do not contribute to RSV-associated ABX


###################################################################################################
#### 3. Calculate the number of averted abx **defined daily doses** (per 1000 person years)    ####
#### attributable to RSV for each age group and each vaccine strategy                          ####
###################################################################################################

###### a) Product of (for each of the n samples and for each of the age groups): 
######## 1. proportion of GP visits averted
######## 2. number of GP visits resulting in ABX course per 100,000 people
######## 3. DDD per course (by age groups) = 7
######## 4. 1,000 / 100,000 (convert from per 100,000 to per 1000)

# ANALYSIS A (Taylor et al. Rx rates), B (Fitzpatrick et al. Rx rate), C (Fitzpatrick et al. Rx rate + no adult Rx)
averted_ddd_analysisA <- averted_ddd_per1000py(gppropaverted, rx_rates_A, popsize)
averted_ddd_analysisB <- averted_ddd_per1000py(gppropaverted, rx_rates_B, popsize)
averted_ddd_analysisC <- averted_ddd_per1000py(gppropaverted, rx_rates_C, popsize)

averted_ddd_0to100yr_per10000py_summary_A <- averted_ddd_analysisA$averted_ddd_0to100yr_per1000py %>%
  pivot_longer(everything(), values_to = "averted_ddd_per1000py", names_to = "intervention") %>%
  group_by(intervention) %>%
  summarise(averted_courses_per10000py_mean = mean(averted_ddd_per1000py * 10 / 7),
            averted_courses_per10000py_loCI = my_quantile(averted_ddd_per1000py * 10 / 7, 0.025),
            averted_courses_per10000py_hiCI = my_quantile(averted_ddd_per1000py * 10 / 7, 0.975)) %>%
            arrange(desc(averted_courses_per10000py_mean))

averted_ddd_0to100yr_per10000py_summary_B <- averted_ddd_analysisB$averted_ddd_0to100yr_per1000py %>%
  pivot_longer(everything(), values_to = "averted_ddd_per1000py", names_to = "intervention") %>%
  group_by(intervention) %>%
  summarise(averted_courses_per10000py_mean = mean(averted_ddd_per1000py * 10 / 7),
            averted_courses_per10000py_loCI = my_quantile(averted_ddd_per1000py * 10 / 7, 0.025),
            averted_courses_per10000py_hiCI = my_quantile(averted_ddd_per1000py * 10 / 7, 0.975)) %>%
  arrange(desc(averted_courses_per10000py_mean))

averted_ddd_0to100yr_per10000py_summary_C <- averted_ddd_analysisC$averted_ddd_0to100yr_per1000py %>%
  pivot_longer(everything(), values_to = "averted_ddd_per1000py", names_to = "intervention") %>%
  group_by(intervention) %>%
  summarise(averted_courses_per10000py_mean = mean(averted_ddd_per1000py * 10 / 7),
            averted_courses_per10000py_loCI = my_quantile(averted_ddd_per1000py * 10 / 7, 0.025),
            averted_courses_per10000py_hiCI = my_quantile(averted_ddd_per1000py * 10 / 7, 0.975)) %>%
  arrange(desc(averted_courses_per10000py_mean))

# save plots of averted courses per 10,000py for each intervention
pa <- abx_reduction_plot(averted_ddd_analysisA, "A")
pb <- abx_reduction_plot(averted_ddd_analysisB, "B")
pc <- abx_reduction_plot(averted_ddd_analysisC, "C")

# print tables of averted courses per 10,000py for each intervention
print(averted_ddd_0to100yr_per10000py_summary_A)
print(averted_ddd_0to100yr_per10000py_summary_B)
print(averted_ddd_0to100yr_per10000py_summary_C)

##### b. Calculate the DDD averted per dose

# create list of total effect by each intervention
### the number of intervention courses
vaccine_courses_mab <- data.frame(inter = c("MAB_VHR_S", "MAB_HR_S","MAB_HR_S+","MAB_ALL_S","MAB_ALL_S+"), 
                              # courses = c(11679, 22907, 22907, 252581, 547818))
                              courses = c(2218, 11679, 22907, 252581, 547818))
vaccine_courses_vax <- data.frame(inter = c("VAC_INF_S", "VAC_INF_A", "VAC_2_4_S","VAC_5_9_S","VAC_5_14_S"),
                              courses =c(251162, 617724, 917008,  2046820,  4093640))
vaccine_courses_mat <- data.frame(inter = c("MAT_A","MAT_S"),
                              courses =  c(165257, 406442))
vaccine_courses <- bind_rows(vaccine_courses_mab, vaccine_courses_vax, vaccine_courses_mat)

A_ddd_0to100_per1000py <- as.list(averted_ddd_analysisA$averted_ddd_0to100yr_per1000py) 
B_ddd_0to100_per1000py <- as.list(averted_ddd_analysisB$averted_ddd_0to100yr_per1000py) 
C_ddd_0to100_per1000py <- as.list(averted_ddd_analysisC$averted_ddd_0to100yr_per1000py) 


averted_ddd_analysisA_per1000doses <- as.data.frame(matrix(unlist(purrr::map(
                              .x = new_interventions,
                               .f = ~(1000 * 1000 * A_ddd_0to100_per1000py[[.x]]) / as.numeric(vaccine_courses[vaccine_courses$inter==.x,"courses"]))
                                      ), ncol = length(new_interventions)))
averted_ddd_analysisB_per1000doses <- as.data.frame(matrix(unlist(purrr::map(
                              .x = new_interventions,
                              .f = ~(1000 * 1000 * B_ddd_0to100_per1000py[[.x]]) / as.numeric(vaccine_courses[vaccine_courses$inter==.x,"courses"]))
                            ), ncol = length(new_interventions)))

averted_ddd_analysisC_per1000doses <- as.data.frame(matrix(unlist(purrr::map(
                              .x = new_interventions,
                              .f = ~(1000 * 1000 * C_ddd_0to100_per1000py[[.x]]) / as.numeric(vaccine_courses[vaccine_courses$inter==.x,"courses"]))
                            ), ncol = length(new_interventions)))

names(averted_ddd_analysisA_per1000doses) <- new_interventions
names(averted_ddd_analysisB_per1000doses) <- new_interventions
names(averted_ddd_analysisC_per1000doses) <- new_interventions


# output efficiencies by converting to courses per person year from DDD per person year
averted_courses_perdose_summary_A <- averted_ddd_analysisA_per1000doses %>%
  pivot_longer(everything(), values_to = "averted_ddd_per1py_per1000doses", names_to = "intervention") %>%
  group_by(intervention) %>%
  summarise(averted_courses_per1py_per1000doses_mean = mean(averted_ddd_per1py_per1000doses / 7),
            averted_courses_per1py_per1000doses_loCI = my_quantile(averted_ddd_per1py_per1000doses / 7, 0.025),
            averted_courses_per1py_per1000doses_hiCI = my_quantile(averted_ddd_per1py_per1000doses / 7, 0.975)) %>%
  arrange(desc(averted_courses_per1py_per1000doses_mean))

# print averted courses per person year per 1,000 vaccine/mAb courses
print(averted_courses_perdose_summary_A)

# save plots of averted DDD per 1000py per vaccine course for each intervention
qa <- abx_reduction_plot_perdose(averted_ddd_analysisA_per1000doses, "A")
qb <- abx_reduction_plot_perdose(averted_ddd_analysisB_per1000doses, "B")
qc <- abx_reduction_plot_perdose(averted_ddd_analysisC_per1000doses, "C")

plot_and_save(list(pa, qa), "A")
plot_and_save(list(pb, qb), "B")
plot_and_save(list(pc, qc), "C")

######################################################################
#### 4. CALCULATE % REDUCTION AS FRACTION OF TOTAL ABX USE        ####
#### using England and Wales Rx data and E&W population sizes     ####
######################################################################

usage_2018 <- load.usage() %>% filter(year == 2018) %>% select(V1) %>% sum()

fractionaverted_oftotalRx_A <- (averted_ddd_analysisA$averted_ddd_0to100yr_per1000py / 7) *
                                          (ew_popsize0to100yr_2018 / 1000) / usage_2018

fractionaverted_oftotalRx_B <- (averted_ddd_analysisB$averted_ddd_0to100yr_per1000py / 7) *
                                          (ew_popsize0to100yr_2018 / 1000) / usage_2018

fractionaverted_oftotalRx_C <- (averted_ddd_analysisC$averted_ddd_0to100yr_per1000py / 7) *
  (ew_popsize0to100yr_2018 / 1000) / usage_2018

###### output table of fraction of ABX that is averted due to RSV vaccine in descending order

fractionaverted_oftotalRx_A %>%
            pivot_longer(everything(), values_to = "ProportionTotalRxAverted", names_to = "intervention") %>%
              group_by(intervention) %>%
              summarise(ProportionTotalRxAverted_mean = mean(ProportionTotalRxAverted),
                        ProportionTotalRxAverted_loCI = my_quantile(ProportionTotalRxAverted, 0.025),
                        ProportionTotalRxAverted_hiCI = my_quantile(ProportionTotalRxAverted, 0.975)) %>%
              arrange(desc(ProportionTotalRxAverted_mean))

fractionaverted_oftotalRx_B %>%
              pivot_longer(everything(), values_to = "ProportionTotalRxAverted", names_to = "intervention") %>%
              group_by(intervention) %>%
              summarise(ProportionTotalRxAverted_mean = mean(ProportionTotalRxAverted),
                        ProportionTotalRxAverted_loCI = my_quantile(ProportionTotalRxAverted, 0.025),
                        ProportionTotalRxAverted_hiCI = my_quantile(ProportionTotalRxAverted, 0.975)) %>%
              arrange(desc(ProportionTotalRxAverted_mean))

fractionaverted_oftotalRx_C %>%
  pivot_longer(everything(), values_to = "ProportionTotalRxAverted", names_to = "intervention") %>%
  group_by(intervention) %>%
  summarise(ProportionTotalRxAverted_mean = mean(ProportionTotalRxAverted),
            ProportionTotalRxAverted_loCI = my_quantile(ProportionTotalRxAverted, 0.025),
            ProportionTotalRxAverted_hiCI = my_quantile(ProportionTotalRxAverted, 0.975)) %>%
  arrange(desc(ProportionTotalRxAverted_mean))
                    
##########################################################################################
#### 5. Calculate the resistance outcomes for averted abx DDD (per 1000 person years) ####
####    attributable to RSV for each age group and each vaccine strategy              ####
##########################################################################################

# Cost calculations
# 1415 in USD 2016 (Shrestha et al.) / 101.38% (Bureau of Labour Statistics CPI) -> USD 2014
costperinfection_2014usd <- 1415 / 1.0138
# = 1395.739 USD 2014 x 79 (GBR) / 130 (USD) health care PPP 2014 (https://www.oecd.org/health/health-systems/
# International-Comparisons-of-Health-Prices-and-Volumes-New-Findings.pdf) x 1.647701 GBP
costperinfection_2014gbp <- costperinfection_2014usd * 79 / 130 * 1 / 1.647701
# = 514.7656 GBP (2014) x 1.14 -> 2020 GBP using https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator
costperinfection_2020gbp <- costperinfection_2014gbp * 1.14
# = 586.8328 GBP 2020

# convert "averted DDD per 1000py" to "gain in XXX per 1000 DDD averted"
res_outcomesA <- calculate_resistanceoutcomes(averted_ddd_analysisA, "A", new_interventions, engpopsize0to100yr, vaccine_courses, costperinfection_2020gbp)
res_outcomesB <- calculate_resistanceoutcomes(averted_ddd_analysisB, "B", new_interventions, engpopsize0to100yr, vaccine_courses, costperinfection_2020gbp)
res_outcomesC <- calculate_resistanceoutcomes(averted_ddd_analysisC, "C", new_interventions, engpopsize0to100yr, vaccine_courses, costperinfection_2020gbp)
