##################################
####   ESTIMATION FUNCTIONS   ####
##################################

# note: si_rate changed from 0.543
# new rate gives average of 7.47 months, or 0.62 years, which is correct


find_rts_exp = function(dat, si_shape = 1.47, si_rate = 0.1968, k = 200,
                        file_name = 'overall_us_exp'){
  # estimate reproductive numbers for overall us incorporating extrapulmonary cases
  
  dates = seq(1, nrow(dat))
  
  # get discrete gamma distribution density
  ps.si = rep(0,k)
  for(i in 1:k){
    ps.si[i] = pgamma(i,shape=si_shape, rate=si_rate) - pgamma(i-1,shape=si_shape,rate=si_rate)
  }
  ps.si = ps.si/sum(ps.si) # this ensures that the ps add to one and make a true distribution
  
  # get weights from serial distribution
  si.weights = matrix(nrow = length(dates), ncol = length(dates))
  for(t1 in 2:length(dates)){
    for(t2 in 1:length(dates)){
      dateDiff = dates[t1] - dates[t2]
      p.si = ifelse(dateDiff>0, ps.si[dateDiff], 0)
      si.weights[t1,t2] = p.si
    }
  }
  si.weights[1,] = 0
  si.weights[is.na(si.weights)] = 0
  
  reweighed.row.sum = si.weights%*%(dat$cases - dat$extra)
  reweighed.prob = si.weights/as.vector(reweighed.row.sum)
  reweighed.prob[is.na(reweighed.prob)] = 0
  
  rts = colSums(reweighed.prob*dat$cases)
  
  to_save <- data.frame(rts, dates)
  colnames(to_save) <- c('rts', 'date')
  
  #write.csv(to_save, file = paste('output/', file_name, '.csv', sep=''), row.names = FALSE)
  return(to_save)
}



find_rts_import = function(data, si_shape = 1.47, si_rate = 0.1968, k = 200, file_name = 'overall_us_import'){
  
  # find reproductive numbers for overall us while incorporating importation (and extrapulmonary)
  dates = seq(1, nrow(data))
  
  # get discrete gamma distribution density
  ps.si = rep(0,k)
  for(i in 1:k){
    ps.si[i] = pgamma(i,shape=si_shape, rate=si_rate) - pgamma(i-1,shape=si_shape,rate=si_rate)
  }
  ps.si = ps.si/sum(ps.si) # this ensures that the ps add to one and make a true distribution
  
  # get weights from serial distribution
  si.weights = matrix(nrow = length(dates), ncol = length(dates))
  for(t1 in 2:length(dates)){
    for(t2 in 1:length(dates)){
      dateDiff = dates[t1] - dates[t2]
      p.si = ifelse(dateDiff>0, ps.si[dateDiff], 0)
      si.weights[t1,t2] = p.si
    }
  }
  si.weights[1,] = 0
  si.weights[is.na(si.weights)] = 0
  
  # weigh by the number of cases
  reweighed.row.sum = si.weights%*%(data$cases - data$extra)
  reweighed.prob = si.weights/as.vector(reweighed.row.sum)
  reweighed.prob[is.na(reweighed.prob)] <- 0
  
  not_imports = data$cases - data$imports
  
  rts <- colSums(reweighed.prob*not_imports)
  
  to_save <- data.frame(rts, data$date)
  colnames(to_save) <- c('rts', 'date')
  
  #write.csv(to_save, file = paste('output/', file_name, '.csv', sep=''),
  #          row.names = FALSE)
  return(to_save)
}

##### estimate Rts with heterogeneity #####

find_rts_strata <- function(data, num_strata = 2, int_mat = diag(2), si_shape = 1.47,
                            si_rate = 0.1968, k = 100, file_name = 'default_name'){
  
  dates <- seq(1, nrow(data) / num_strata)
  
  # get discrete gamma distribution density
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
      pgamma(i - 1, shape = si_shape, rate = si_rate)
  }
  ps_si <- ps_si / sum(ps_si)
  
  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0
  
  
  # reweigh these weights by the interaction between strata
  int_weights = si_weights %x% int_mat
  
  # weigh by the number of cases
  reweighed_row_sum <- int_weights %*% (data$cases - data$extra)
  reweighed_prob <- int_weights / as.vector(reweighed_row_sum)
  reweighed_prob[is.na(reweighed_prob)] <- 0
  
  rts_raw <- colSums(reweighed_prob * data$cases)  # calculate reproductive numbers
  
  to_save <- data.frame(rts_raw, data$date)
  colnames(to_save) <- c('rts', 'date')
  
  #write.csv(to_save, file = paste('output/', file_name, '.csv', sep=''), row.names = FALSE)
  return(to_save)
  
}


find_rts_strata_2 = function(data, num_strata = 2, int_mat = diag(2), num_strata_2 = 7, int_mat_2 = diag(7),
                             si_shape = 1.47, si_rate = 0.1968, k = 200, file_name = 'default_name'){
  
  dates <- seq(1, nrow(data) / (num_strata * num_strata_2))
  
  # get discrete gamma distribution density
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
      pgamma(i - 1, shape = si_shape, rate = si_rate)
  }
  ps_si <- ps_si / sum(ps_si)
  
  
  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0
  
  # reweigh these weights by the interaction between strata
  int_weights = si_weights %x% int_mat %x% int_mat_2
  
  # weigh by the number of cases
  reweighed_row_sum <- int_weights %*% (data$cases - data$extra)
  reweighed_prob <- int_weights / as.vector(reweighed_row_sum)
  
  reweighed_prob[is.na(reweighed_prob)] <- 0
  
  rts_raw <- colSums(reweighed_prob * data$cases)  # calculate reproductive numbers
  
  to_save <- data.frame(rts_raw, data$date)
  colnames(to_save) <- c('rts', 'date')
  
  # write.csv(to_save, file = paste('output/', file_name, '.csv', sep=''),
  #           row.names = FALSE)
  
  return(to_save)
  
}


# find rts with one category (eg. region) and importation
find_rts_complete = function(data, num_strata = 12, int_mat = diag(12), si_shape = 1.47,
                             si_rate = 0.1968, k = 200, file_name = 'default_name'){
  
  dates <- seq(1, nrow(data) / num_strata)
  
  
  # get discrete gamma distribution density
  ps.si = rep(0,k)
  for(i in 1:k){
    ps.si[i] = pgamma(i,shape=si_shape, rate=si_rate) - pgamma(i-1,shape=si_shape,rate=si_rate)
  }
  ps.si = ps.si/sum(ps.si) # this ensures that the ps add to one and make a true distribution
  
  
  # get weights from serial distribution
  si.weights = matrix(nrow = length(dates), ncol = length(dates))
  for(t1 in 2:length(dates)){
    for(t2 in 1:length(dates)){
      dateDiff = dates[t1] - dates[t2]
      p.si = ifelse(dateDiff>0, ps.si[dateDiff], 0)
      si.weights[t1,t2] = p.si
    }
  }
  si.weights[1,] = 0
  si.weights[is.na(si.weights)] = 0
  
  # reweigh these weights by the interaction between regions
  int.weights = si.weights %x% int_mat
  
  # weigh by the number of cases
  reweighed.row.sum = int.weights %*% (data$cases - data$extra)
  reweighed.prob = int.weights/as.vector(reweighed.row.sum)
  reweighed.prob[is.na(reweighed.prob)] = 0
  
  not_imports = data$cases - data$imports
  
  rts = colSums(reweighed.prob * not_imports)
  
  to_save <- data.frame(rts, data$date)
  colnames(to_save) <- c('rts', 'date')
  
  # write.csv(to_save, file = paste('output/', file_name, '.csv', sep=''),
  #           row.names = FALSE)
  
  return(to_save)
}


find_rts_complete_2 = function(data, num_strata = 2, int_mat = diag(2), num_strata_2 = 7, int_mat_2 = diag(7),
                               si_shape = 1.47, si_rate = 0.1968, k = 200, file_name = 'default_name'){
  
  dates <- seq(1, nrow(data) / (num_strata * num_strata_2))
  
  # get discrete gamma distribution density
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
      pgamma(i - 1, shape = si_shape, rate = si_rate)
  }
  ps_si <- ps_si / sum(ps_si)
  
  
  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0
  
  # reweigh these weights by the interaction between strata
  int_weights = si_weights %x% int_mat %x% int_mat_2
  
  # weigh by the number of cases
  reweighed_row_sum <- int_weights %*% (data$cases - data$extra)
  reweighed_prob <- int_weights / as.vector(reweighed_row_sum)
  
  reweighed_prob[is.na(reweighed_prob)] <- 0
  
  not_imports = data$cases - data$imports
  
  rts = colSums(reweighed_prob * not_imports)
  
  to_save <- data.frame(rts, data$date)
  colnames(to_save) <- c('rts', 'date')
  
  # write.csv(to_save, file = paste('output/', file_name, '.csv', sep=''),
  #           row.names = FALSE)
  
  return(to_save)
  
}



##########################
##      DATA SETUP      ##
##########################

# set up data for overall US
setup_import = function(data){
  
  dat_g = data %>% group_by(date)
  
  case_sum = dat_g %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total)
  extra_sum = dat_g %>% filter(extra == 'EXTRAPULM ONLY') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total) 
  import_sum = dat_g %>% filter(rt == '0') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total) 
  
  dat = data.frame(case_sum$date, case_sum$cases_total, extra_sum$cases_total, import_sum$cases_total)
  colnames(dat) = c('date', 'cases', 'extra', 'imports')
  
  return(dat)
}

# set up data for overall us with US born categories
setup_usborn = function(data){
  
  dat_g = data %>% group_by(date, usborn)
  
  case_sum = dat_g %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, usborn)
  extra_sum = dat_g %>% filter(extra == 'EXTRAPULM ONLY') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, usborn) 
  import_sum = dat_g %>% filter(rt == '0') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, usborn) 
  
  dat = data.frame(case_sum$date, case_sum$usborn, case_sum$cases_total, extra_sum$cases_total, import_sum$cases_total)
  colnames(dat) = c('date', 'usborn', 'cases', 'extra', 'imports')
  
  return(dat)
}

# overall us with age
setup_age = function(data){
  
  dat_g = data %>% group_by(date, age)
  
  case_sum = dat_g %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, age)
  extra_sum = dat_g %>% filter(extra == 'EXTRAPULM ONLY') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, age) 
  import_sum = dat_g %>% filter(rt == '0') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, age) 
  
  dat = data.frame(case_sum$date, case_sum$age, case_sum$cases_total, extra_sum$cases_total, import_sum$cases_total)
  colnames(dat) = c('date', 'age', 'cases', 'extra', 'imports')
  
  return(dat)
}


# overall us with age and US born categories
setup_age_usborn = function(data){
  
  dat_g = data %>% group_by(date, age, usborn)
  
  case_sum = dat_g %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, age, usborn)
  extra_sum = dat_g %>% filter(extra == 'EXTRAPULM ONLY') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, age, usborn) 
  import_sum = dat_g %>% filter(rt == '0') %>% summarize(cases_total = sum(cases)) %>% ungroup() %>% select(date, cases_total, age, usborn) 
  
  dat = data.frame(case_sum$date, case_sum$age, case_sum$usborn, case_sum$cases_total, extra_sum$cases_total, import_sum$cases_total)
  colnames(dat) = c('date', 'age', 'usborn', 'cases', 'extra', 'imports')
  
  return(dat)
}



# bootstrapping
bootstrap_samples = function(data, n = 1000){
  # n = number of iterations
  # data = processed data
  
  n = 1000  # number of iterations
  all_iter = list()
  for (j in 1:n){
    
    # block size
    b = 4
    
    date_list = unique(data$date)
    # get indices where to start blocks
    blocks = sample((length(date_list)-b), size = floor(length(date_list)/b + 1), replace = TRUE)
    blocks_sort = sort(blocks)
    
    boot = data.frame()
    
    # create bootstrap samples
    for (i in 1:(length(date_list)/b + 1)){  # change to length(blocks)
      date_block = date_list[blocks_sort[i]:(blocks_sort[i]+b-1)]
      full_block =  data[which(data$date %in% date_block==T),]
      boot = rbind(boot, full_block)
    }
    
    # cut it to the length of the original data
    cut_boot = boot[1:(nrow(data)), ]
    
    all_iter[[j]] = cut_boot
  }
  
  return(all_iter)
}


bootstrap_estimates = function(samples, est_option){
  # samples is output from bootstrap_samples
  # estimate_option - specify which estimation function to use
  n = length(samples)
  rt_df = data.frame(matrix(ncol=n, nrow=nrow(samples[[1]])))
  
  if (est_option == 'usborn_age_import'){
    for(i in 1:1000){
      rt = find_rts_complete_2(samples[[i]], num_strata = 3, int_mat = diag(3),
                               num_strata_2 = 7, int_mat_2 = diag(7))
      rt_df[i] = rt$rts
    }
  }
  
  else if (est_option == 'usborn_age_import_2'){
    for(i in 1:1000){
      rt = find_rts_complete_2(samples[[i]], num_strata = 3, int_mat = diag(3),
                               num_strata_2 = 7, int_mat_2 = age_matrix)
      rt_df[i] = rt$rts
    }
  }
  
  else if (est_option == 'all'){
    for(i in 1:1000){
      rt = find_rts_strata(samples[[i]], num_strata = 1, int_mat = diag(1))
      rt_df[i] = rt$rts
    }
  }
  
  else if (est_option == 'all_import'){
    for(i in 1:1000){
      rt = find_rts_complete(samples[[i]], num_strata = 1, int_mat = diag(1))
      rt_df[i] = rt$rts
    }
  }
  
  else if (est_option == 'all_usborn_import'){
    for(i in 1:1000){
      rt = find_rts_complete(samples[[i]], num_strata = 3, int_mat = diag(3))
      rt_df[i] = rt$rts
    }
  }
  
  else if (est_option == 'all_age_import'){
    for(i in 1:1000){
      rt = find_rts_complete(samples[[i]], num_strata = 7, int_mat = diag(7))
      rt_df[i] = rt$rts
    }
  }
  
  else if (est_option == 'all_age_import_2'){
    for(i in 1:1000){
      rt = find_rts_complete(samples[[i]], num_strata = 7, int_mat = age_matrix)
      rt_df[i] = rt$rts
    }
  }
  
  return(rt_df)
}


get_ci = function(estimates){
  # estimates is dataframe of bootstrap estimates
  ci = data.frame(matrix(ncol=7, nrow=nrow(estimates)))
  for(i in 1:nrow(estimates)){
    rowi = as.numeric(estimates[i,])
    
    # ci estimates
    lower = quantile(rowi, 0.025, na.rm = TRUE)
    upper = quantile(rowi, 0.975, na.rm = TRUE)
    ci[i,1] = lower
    ci[i,2] = upper
    
    # for boxplots
    mini = min(rowi, na.rm = TRUE)
    maxi = max(rowi, na.rm = TRUE)
    medi = median(rowi, na.rm = TRUE)
    q25 = quantile(rowi, 0.25, na.rm = TRUE)
    q75 = quantile(rowi, 0.75, na.rm = TRUE)
    
    ci[i,3] = mini
    ci[i,4] = maxi
    ci[i,5] = medi
    ci[i,6] = q25
    ci[i,7] = q75
    
  }
  colnames(ci) = c('lower','upper','min','max','median','q25','q75')
  
  return(ci)
}


avg_est_ci = function(est, est_boot, time = 'all', categories = 'none'){
  # est: estimated reproductive number dataframe
  # est_boot: bootstrap estimate dataframe
  # time: 'all' for entire time span (2011-2015, cuts out first year and last two years),
  # 'year' to get estimate by year (again, 2011-2015)
  # categories: categories to stratify by
  
  boot0 = data.frame(est, est_boot)
  boot0$year = substr(boot0$date, 1, 4)
  
  # cut off first year and last two years
  boot = boot0 %>%
    filter(year %in% 2011:2017)
  
  
  if(time == 'all'){
    
    rt_est = mean(boot$rts, na.rm = TRUE)
    
    if(categories == 'none'){
      boot_grouped = boot %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE)
      rt_est = boot %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'usborn'){
      boot_grouped = boot %>%
        group_by(usborn) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE)
      rt_est = boot %>%
        group_by(usborn) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'age'){
      boot_grouped = boot %>%
        group_by(age) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE)
      rt_est = boot %>%
        group_by(age) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'region'){
      boot_grouped = boot %>%
        group_by(region) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE)
      rt_est = boot %>%
        group_by(region) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'age, region'){
      boot_grouped = boot %>%
        group_by(age, region) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE)
      rt_est = boot %>%
        group_by(age, region) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
  }
  
  if(time == 'year'){
    
    if(categories == 'none'){
      boot_grouped = boot %>%
        group_by(year) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE) %>%
        ungroup()
      rt_est = boot %>%
        group_by(year) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'usborn'){
      boot_grouped = boot %>%
        group_by(year, usborn) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE) %>%
        ungroup()
      rt_est = boot %>%
        group_by(year, usborn) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'age'){
      boot_grouped = boot %>%
        group_by(year, age) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE) %>%
        ungroup()
      rt_est = boot %>%
        group_by(year, age) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'region'){
      boot_grouped = boot %>%
        group_by(year, region) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE) %>%
        ungroup()
      rt_est = boot %>%
        group_by(year, region) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
    if(categories == 'age, region'){
      boot_grouped = boot %>%
        group_by(year, age, region) %>%
        summarise_at(vars(X1:X1000), mean, na.rm = TRUE) %>%
        ungroup()
      rt_est = boot %>%
        group_by(year, age, region) %>%
        summarise(rt_mean = mean(rts, na.rm = TRUE))
    }
    
  }
  
  for_ci = boot_grouped %>%
    select(X1:X1000)
  ci = get_ci(for_ci)
  
  labels = boot_grouped %>%
    select(-starts_with('X'))
  
  ci_labeled = data.frame(rt_est$rt_mean, labels, ci)
  colnames(ci_labeled)[1] = 'rt_est' 
  
  return(ci_labeled)
}










