# FUNCTIONS
### TAKEN / ADAPTED FROM CHAE et al. (2018) Emerging Infectious Diseases
# function to calculate quantiles
my_quantile <- function(x, probs) {
  dplyr::tibble(x = quantile(x, probs, na.rm = TRUE), probs = probs)$x
}

# Draw from triangular distribution specified according to mode and HDI (l,r)
# containing proportion p of density, clamped to range [minimum, maximum]
rtriang <- function(n, mode, l, r, p = 0.95, minimum = -Inf, maximum = Inf) {
  # bounds of triangular distribution: a = left vertex, b = right vertex, c = peak
  h <- sqrt(1 - p);
  a <- (h * mode - l) / (h - 1);
  b <- (r - h * mode) / (1 - h);
  c <- mode;
  
  # choose variates from triangular distribution
  u <- runif(n);
  x <- ifelse(u < (c - a) / (b - a), a + sqrt(u * (b - a) * (c - a)), b - sqrt((1 - u) * (b - a) * (b - c)));
  return (pmin(maximum, pmax(minimum, x)));
}

# What does this do? Calculate the cases averted per daily dose per year?
averted_ddd_per1000py <- function(gp_propaverted, rx_rates, popsize) {
  
  averted_ddd_0to5mo_per1000py <- 0.01 * 7 * gp_propaverted$GPvisits_proportionaverted_0to5mo * rx_rates$Rx_0to5mo_per100000py
  averted_ddd_6to23mo_per1000py <- 0.01 * 7 * gp_propaverted$GPvisits_proportionaverted_6to23mo * rx_rates$Rx_6to23mo_per100000py
  averted_ddd_2to4yr_per1000py <- 0.01 * 7 * gp_propaverted$GPvisits_proportionaverted_2to4yr * rx_rates$Rx_2to4y_per100000py
  averted_ddd_5to17yr_per1000py <- 0.01 * 7 * gp_propaverted$GPvisits_proportionaverted_5to17yr * rx_rates$Rx_5to17y_per100000py
  averted_ddd_18to100yr_per1000py <- 0.01 * 7 * gp_propaverted$GPvisits_proportionaverted_18to100yr * rx_rates$Rx_18to100y_per100000py
  
  # combine all age groups to find averted ddd per 1000py across all groups
  averted_ddd_0to100yr_per1000py <- (1 / popsize$engpopsize0to100yr) *
    (averted_ddd_0to5mo_per1000py * popsize$engpopsize0to5mo
     + averted_ddd_6to23mo_per1000py * popsize$engpopsize6to23mo
     + averted_ddd_2to4yr_per1000py * popsize$engpopsize2to4yr
     + averted_ddd_5to17yr_per1000py * popsize$engpopsize5to17yr
     + averted_ddd_18to100yr_per1000py * popsize$engpopsize18to100yr)

  return(list(averted_ddd_0to5mo_per1000py = averted_ddd_0to5mo_per1000py,
              averted_ddd_6to23mo_per1000py = averted_ddd_6to23mo_per1000py,
              averted_ddd_2to4yr_per1000py = averted_ddd_2to4yr_per1000py,
              averted_ddd_5to17yr_per1000py = averted_ddd_5to17yr_per1000py,
              averted_ddd_18to100yr_per1000py = averted_ddd_18to100yr_per1000py,
              averted_ddd_0to100yr_per1000py = averted_ddd_0to100yr_per1000py))
  
}

# What does this do? 
abx_reduction_plot <- function(averted_ddd_per1000py, analysis){
  abxreduction <- list("0-5mo" = averted_ddd_per1000py$averted_ddd_0to5mo_per1000py, 
                       "6-23mo" = averted_ddd_per1000py$averted_ddd_6to23mo_per1000py,
                       "2-4yr" = averted_ddd_per1000py$averted_ddd_2to4yr_per1000py,
                       "5-17yr" = averted_ddd_per1000py$averted_ddd_5to17yr_per1000py,
                       "18-100yr" = averted_ddd_per1000py$averted_ddd_18to100yr_per1000py,
                       "0-100yr" = averted_ddd_per1000py$averted_ddd_0to100yr_per1000py) %>%
    dplyr::bind_rows(.id = "age_group") %>%
    tidyr::pivot_longer(-age_group, values_to = "averted_ddd_per1000py", names_to = "intervention") %>%
    dplyr::arrange(averted_ddd_per1000py) %>%
    dplyr::mutate(age_group = factor(age_group, levels=c("0-5mo", "6-23mo", "2-4yr", "5-17yr", "18-100yr", "0-100yr"))) %>%
    dplyr::mutate(intervention = factor(intervention, levels=c("MAB_VHR_S", "MAB_HR_S", "MAB_HR_S+", "MAB_ALL_S", "MAB_ALL_S+", "MAT_S", "MAT_A",     
                                                            "VAC_INF_S","VAC_INF_A" , "VAC_2_4_S", "VAC_5_9_S", "VAC_5_14_S")))
  
  
  p <- ggplot2::ggplot(data = abxreduction) + 
    # ggplot2::geom_violin(aes(x = intervention, y = averted_ddd_per1000py, color = age_group)) +
    ggplot2::geom_boxplot(aes(x = intervention, y = averted_ddd_per1000py, color = age_group)) +
    theme_bw() +
    facet_wrap(~age_group, nrow = 3, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    theme(legend.position = "none") +
    ylab("Averted DDD per 1000 person years") +
    xlab("")
  
  return(p)
}

# What does this do?
abx_reduction_plot_perdose <- function(averted_ddd_per1000py_perdose, analysis){
  
  abxreduction <- averted_ddd_per1000py_perdose %>%
    tidyr::pivot_longer(everything(), values_to = "averted_ddd_per1000py", names_to = "intervention")  %>%
    dplyr::mutate(intervention = factor(intervention, levels=c("MAB_VHR_S", "MAB_HR_S", "MAB_HR_S+", "MAB_ALL_S", "MAB_ALL_S+", "MAT_S", "MAT_A",     
                                                            "VAC_INF_S","VAC_INF_A" , "VAC_2_4_S", "VAC_5_9_S", "VAC_5_14_S")))
  
  p <- ggplot2::ggplot(data = abxreduction) + 
    ggplot2::geom_boxplot(aes(x = intervention, y = averted_ddd_per1000py)) +
    # facet_wrap(~age_group, nrow = 3, scales = "free_y") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=  1)) +
    theme(legend.position = "none") +
    ylab("Averted DDD per year \n per 1000 courses") +
    xlab("")
  return(p)
}

plot_and_save <- function(plot_list, analysis){

  P <- grid.arrange(plot_list[[1]], plot_list[[2]],
                    nrow = 2,
                    # layout_matrix = rbind(c(1,1), c(2,3)),
                    # widths = c(2.7, 2.7),
                    heights = c(8.5 * (2 / 3), 8.5 * (1 / 3)))
  todaysdate <- format(Sys.Date(), "%Y%m%d")
        
        ggplot2::ggsave(
          here::here("figs_averted_ddd", paste("averted_ddd_", analysis, "_", todaysdate, ".pdf", sep="")),
          plot = P,
          width = 11,
          height = 8.5,
          units = "in",
          dpi = 300) 
        
        ggplot2::ggsave(
          here::here("figs_averted_ddd", paste("averted_ddd_", analysis, "_", todaysdate, ".png", sep="")),
          plot = P,
          device = png(),
          width = 11,
          height = 8.5,
          units = "in",
          dpi = 300) 
        
    P
}


predict <- function(data, countries, strains, outcomes, predictors, projections, pop.factor, intervention){
  # predict: predict the impact of reducing prescribing on outcomes. 
  # data is a data.table from the function load.data;
  # countries is the vector of countries;
  # strains is a vector of strains to predict the impact upon;
  # outcomes is a vector of outcomes (e.g. DALY.mid);
  # predictors is a vector of predictors (e.g. J01_AC_2015);
  # projections is a vector of the projected change in each predictor;
  
  impacts = NULL;
  for (outcome in outcomes) {
    for (str in strains) {
      # Build model
      subdata = data[strain == str, c(..outcome, ..predictors)];
      model = lm(paste(outcome, "~ ."), data = subdata);
      
      # Make predictions
      co = countries
      # for (co in countries) {
      # change = 0;
      # for (p in 1:length(predictors)) {
      # change = change + projections[p] * model$coeff[p + 1] * 
      #   data[country == co & strain == str, norm];
      
      # calculate the change for each of the projections
      # change = projections * model$coeff[2] * 
      #   data[country == co & strain == str, norm];
      
      change = projections * model$coeff[2] * 
        data[country == co & strain == str, norm];
      
      
      # }
      impact = data.table( #country = co, 
        strain = str, 
        outcome = outcome, 
        # predictors = paste(predictors, collapse = "+"),
        # orig = data[country == co & strain == str][[outcome]] * data[country == co & strain == str, norm] * pop.factor,
        change = change * pop.factor,
        iter = 1:n,
        inter = intervention)
      # final = data[country == co & strain == str, ..outcome][[outcome]] * 
      # data[country == co & strain == str, norm] * pop.factor + change * pop.factor);
      
      impacts = rbind(impacts, impact);
      # }
    }
  }
  return (impacts)
}

calculate_resistanceoutcomes <- function(averted_ddd, analysis, new_interventions, popsize, vaccine_courses, drinfection_cost){
  
  reg <- load.data(2015, "AC", "J01", "BSI", max);
  impacts <- NULL
  impact <- NULL
  # age.dist = fread("laiv_amr_ew/ONS_EW_agedist.csv")
  # pop = age.dist$population_2015 # use the 
  for (intervention in new_interventions){
    # grab the total gain in outcomes
    
    effect_size <- -averted_ddd$averted_ddd_0to100yr_per1000py[[intervention]] / 365
    impact <- predict(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                      "J01_AC_2015", effect_size, (popsize / 65128861), intervention)
    impacts <- rbind(impacts, impact)
    
  }
  
  impact_all <- impacts[, .(change = -sum(change)), by = list(outcome, iter, inter)]
  
  impact_all$outcome <- impact_all$outcome %>%
                    recode_factor(DALY.mid = "DALYs", cases.mid = "Cases", deaths.mid = "Deaths")
  impact_all <- impact_all %>%
    dplyr::mutate(inter = factor(inter, levels=c("MAB_VHR_S", "MAB_HR_S", "MAB_HR_S+", "MAB_ALL_S", "MAB_ALL_S+", "MAT_S", "MAT_A",     
                                                               "VAC_INF_S","VAC_INF_A" , "VAC_2_4_S", "VAC_5_9_S", "VAC_5_14_S"))) %>%
    rowwise() %>%
    dplyr::mutate(change_perdose = 1000 * change / vaccine_courses[vaccine_courses$inter == inter,"courses"] ) %>%
    dplyr::mutate(cost = ifelse(outcome == "Cases", drinfection_cost * change, NA))
  
  print( impact_all %>%  group_by(outcome, inter) %>%
    summarise(mean = mean(change),
              loCI = my_quantile(change, 0.025) ,
              hiCI = my_quantile(change, 0.975),
              mean_change_per1000doses = mean(change_perdose),
              mean_cost = mean(cost),
              cost_loCI = my_quantile(cost, 0.025) ,
              cost_hiCI = my_quantile(cost, 0.975)) %>%
      arrange(outcome, desc(mean)), n = 36)
  
  p_o <- ggplot2::ggplot(data = impact_all) + 
    ggplot2::geom_boxplot(aes(x = inter, y = change)) +
    facet_wrap(~outcome, nrow = 3, scales = "free_y") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    theme(legend.position = "none") +
    ylab("Annual averted outcomes") +
    xlab("")
  
  todaysdate <- format(Sys.Date(), "%Y%m%d")
  
  ggplot2::ggsave(
    here::here("figs_averted_pop_outcomes", paste("averted_pop_outcomes_", analysis, "_", todaysdate, ".png", sep="")),
    plot = p_o,
    width = 11,
    height = 8.5,
    units = "in",
    dpi = 300) 
  
  ggplot2::ggsave(
    here::here("figs_averted_pop_outcomes", paste("averted_pop_outcomes_", analysis, "_", todaysdate, ".pdf", sep="")),
    plot = p_o,
    device = png(),
    width = 11,
    height = 8.5,
    units = "in",
    dpi = 300) 
  
  return(impact_all)
  
}

# Printable table
PrintAgeStratified = function(agestrat){
  agestrat = agestrat
  
  table1 = paste0("\t", paste0(c(names(agestrat)[1:2], "% RSV-associated GP visits averted across 10y"), collapse = "\t"), "\n")
  prec = 2
  
  for (int in unique(agestrat$intervention)){
    
    # 
        for (r in c(1,5,3,4,2)){
          data = agestrat %>%
            filter(intervention == int) %>%
            filter(age_group == age_group[r])
          # for (co in 3:ncol(agestrat)) {
            row = paste0(int, "\t",
                         data$age_group, 
                         "\t", as.numeric(signif(data[3], prec)),
                         " (", as.numeric(signif(data[4], prec)),
                         " – ", as.numeric(signif(data[5], prec)), ")");
          # }
          table1 = paste0(table1, row, "\n");
         }
    # 
    table1 = paste0(table1, "\n")
   }
        cat(table1)
}

# load.drugs: loads drug consumption for a subset of years, sectors (AC, HC, or ACHC), and ATC codes.
load.drugs = function(years, sectors, atcs){
  # drugs = data.table::fread(paste0(root, "drugs_2015.txt"));
  drugs = fread("laiv_amr_ew/drugs_2015.txt");
  drugs = drugs[year %in% years & sector %in% sectors & atc %in% atcs];
  return (dcast(drugs, country ~ atc + sector + year, value.var = "ddd.per.thousand"))
}

# load burden: loads resistant disease burdens normalised according to norm (none, pop for
# population by country, or BSI for total bloodstream infections by each species). BSI.agg
# should be either mean or max. Note that this loads the number of cases from original
# files, so the calculated incidences are not standardised for different age distributions.
load.burden = function(norm, BSI.agg = max){
  # burden = fread(paste0(root, "burden.txt"));
  burden = fread(here::here("laiv_amr_ew","burden.txt"));
  
  if (norm == "none") {
    # no normalization
    burden$norm = 1;
  } else if (norm == "pop") {
    # normalize disease burden relative to population in thousands
    # pop.cov = fread(paste0(root, "population_and_coverage.txt"));
    pop.cov = fread(here::here("laiv_amr_ew", "population_and_coverage.txt"));
    burden$norm = pop.cov$pop.2015.eurostat[match(burden$country, pop.cov$country)] / 1000;
  } else if (norm == "BSI") {
    # normalise disease burden relative to total number of bloodstream infections
    # caused by the species in question (whether resistant or sensitive)
    # pop.cov = fread(paste0(root, "population_and_coverage.txt"));
    pop.cov = fread(here::here("laiv_amr_ew", "population_and_coverage.txt"));
    # isolates = fread(paste0(root, "isolates_2015.csv"));
    isolates = fread(here::here("laiv_amr_ew", "isolates_2015.csv"));
    
    # Collate total tested isolates and coverage for all strains
    spp = c("Acinetobacter/ColRACI/CRACI/MDRACI",
            "Enterococcus faecalis/VRE",
            "Enterococcus faecium/VRE",
            "Escherichia/ColREC/CREC/3GCREC",
            "Klebsiella/ColRKP/CRKP/3GCRKP",
            "Pseudomonas/ColRPA/CRPA/MDRPA",
            "Staphylococcus/MRSA",
            "Streptococcus/PRSP/PMRSP");
    bsi = NULL;
    for (s in spp) {
      args = strsplit(s, "/")[[1]];
      species = args[1];
      for (strain in args[2:length(args)]) {
        bsi.uncorr = isolates[Population %like% species,
                              .(strain = ..strain, tested = BSI.agg(NumValue)),
                              by = .(country = RegionName)];
        bsi.correction = cbind.data.frame(country = pop.cov$country, coverage = pop.cov[[paste0("coverage.", strain)]]);
        bsi = rbind(bsi, merge(bsi.uncorr, bsi.correction, by = "country"));
      }
    }
    
    # Sum together both Enterococcus species
    bsi = rbind(bsi[strain != "VRE"],
                bsi[strain == "VRE", .(tested = sum(tested)), by = .(country, strain, coverage)][, c(1,2,4,3)]);
    bsi[, norm := tested / (coverage / 100)];
    
    burden = merge(burden, bsi[, .(country, strain, norm)], by = c("country", "strain"))
  } else {
    stop("norm must be one of: none, pop, BSI.");
  }
  
  # perform normalisation
  burden[, 3:11] = burden[, 3:11] / burden$norm;
  
  # Greece's unreported S. pneumoniae burdens are given as 0 in this data set, so remove these
  burden[country == "Greece" & strain %like% "^PM?RSP$", 3:11] = NA;
  
  return (burden)
}

# load.data: loads both drug consumption and burden data.
load.data = function(years, sectors, atcs, norm, BSI.agg = max){
  d = load.drugs(years, sectors, atcs);
  b = load.burden(norm, BSI.agg);
  return (merge(b, d, by = "country"))
}

load.usage = function(){
  usage = fread(here::here("laiv_amr_ew", "england_bnf0501_total.csv"));
  usage[, year := floor(PERIOD/100)]
  usage[, month := PERIOD - year * 100]
  usage[, time := year + (month-1)/12]
  return(usage)
}   
