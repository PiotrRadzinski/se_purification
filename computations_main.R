library(dplyr)
library(stringr)

## data read
data = do.call(rbind, lapply(list.files(path = 'data/'), function(file){
    cbind('sample' = str_sub(file, end = -5), read.csv(paste0('data/', file)))
}))
data[which(data$sample == 'artificial_sea_water'), 1] = 'ASW'
data[which(data$sample == 'NIST_SRM_3149'), 1] = 'NIST'

## delta computation
data[, 'delta'] = 1e3 * (2 * data$spl_r / (data$std1_r + data$std2_r) - 1)

## propagation of uncertainty computations
f_abs = 2 * data$spl_r / (data$std1_r + data$std2_r)
f_spl = data$spl_se / data$spl_r
f_std = sqrt(data$std1_se^2 + data$std2_se^2) / (data$std1_r + data$std2_r)

data$delta_se = 1e3 * abs(f_abs) * sqrt(f_spl^2 + f_std^2)

## adding assumed delta true values
data[which(data$sample == 'ASW'), 'delta_true'] = 0
data[which(data$sample == 'BC210a'), 'delta_true'] = data$delta[which(data$sample == 'BC210a')] %>% mean()
data[which(data$sample == 'MAG1'), 'delta_true'] = 0.21
data[which(data$sample == 'NASS4'), 'delta_true'] = 0
data[which(data$sample == 'NIST'), 'delta_true'] = 0
data[which(data$sample == 'SCo1'), 'delta_true'] = 0.175
data[which(data$sample == 'SELM1'), 'delta_true'] = -0.68
data[which(data$sample == 'SGR1'), 'delta_true'] = 0.21 # 0.32

## Shapiro-Wilk normality tests
data$deviation = data$delta - data$delta_true

sapply(unique(data$sample), function(j){
    data$deviation[which(data$sample == j)] %>% shapiro.test()
})
data$deviation %>% shapiro.test()

## one-sided Student's t-tests for expected value == 0 
sapply(unique(data$sample), function(j){
    data$deviation[which(data$sample == j)] %>% t.test()
})
data$deviation %>% t.test()

## direct approach -- sigma_delta & quantile
computations = function(df){
    sd_ = sd(df$deviation) 
    qnorm_ = qnorm(0.975, sd = sd_) 
    c('SD' = sd_, 'quantile 0.975' = qnorm_) %>% round(4)
}

sapply(unique(data$sample), function(j){
    data[which(data$sample == j), ] %>% computations()
})
data %>% computations()

## statistics computations
d_avg = aggregate(data$delta, list(data$sample), FUN=mean) 
d_avg[, 2] = round(d_avg[, 2], 4)
d_avg # delta_avg

D_avg = aggregate(data$deviation, list(data$sample), FUN=mean) 
D_avg[, 2] = round(D_avg[, 2], 4)
D_avg # Delta_avg
mean(data$deviation) %>% round(4)

aggregate(data$delta, list(data$sample), FUN=length) # n

# ## sigma_delta & quantile for delta_avg in place of delta_true
# means_by_group = (data %>% group_by(sample) %>% mutate(t=mean(delta)))[, 't'] %>% unlist() %>% unname()
# 
# sd(data$delta - means_by_group) # SD
# qnorm(0.975, sd = sd(data$delta - means_by_group)) # 0.975 quantile

## monte carlo simulation
mc_comps = function(df, n = 1e6){
    k = nrow(df)
    sapply(1:n, function(j) sd(df$delta + rnorm(k, 0, df$delta_se) - df$delta_true)) %>% mean()
}

mc_sym = sapply(unique(data$sample), function(j){
    data[which(data$sample == j), ] %>% mc_comps()
})
mc_sym %>% round(4)
# ASW     BC210a       MAG1      NASS4       NIST       SCo1      SELM1       SGR1 
# 0.08481583 0.08110133 0.16206990 0.10347354 0.08732488 0.14696945 0.08159445 0.07984867 
qnorm(0.975, sd = mc_sym) %>% round(4)

data %>% mc_comps() # 0.09634643
qnorm(0.975, sd = 0.09634643)
