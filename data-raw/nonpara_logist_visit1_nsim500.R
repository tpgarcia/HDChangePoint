
load("C:/Users/unkyung/Dropbox/logist/nonpara/logist_visit1/combined_nsim500_logist_visit1.RData")

nonpara_logist_visit1_nsim500<-combi.res
save(nonpara_logist_visit1_nsim500, file="Data/nonpara_logist_visit1_nsim500.RData")
usethis::use_data(nonpara_logist_visit1_nsim500, overwrite = TRUE)


