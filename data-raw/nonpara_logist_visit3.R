load("C:/Users/unkyung/Dropbox/Project_Posdoc_Tanya/HDproject_local/Revision/code_for_R_package/nonpara_logist_visit3.RData")
nonpara_logist_visit3<-results
save(nonpara_logist_visit3, file="Data/nonpara_logist_visit3.RData")
usethis::use_data(nonpara_logist_visit3, overwrite = TRUE)
