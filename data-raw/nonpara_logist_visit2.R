load("C:/Users/unkyung/Dropbox/Project_Posdoc_Tanya/HDproject_local/Revision/code_for_R_package/nonpara_logist_visit2.RData")
nonpara_logist_visit2<-results
save(nonpara_logist_visit2, file="Data/nonpara_logist_visit2.RData")
usethis::use_data(nonpara_logist_visit2, overwrite = TRUE)
