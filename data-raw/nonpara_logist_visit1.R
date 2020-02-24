## code to prepare `DATASET` dataset goes here
#usethis::use_data("DATASET")

load("C:/Users/unkyunglee/Dropbox/Project_Posdoc_Tanya/HDproject_local/Rcode(Revision_Rpackage)/code_for_R_package/HDChangePoint/Data/nonpara_logist_visit1.RData")

nonpara_logist_visit1<-combined_results
save(nonpara_logist_visit1, file="Data/nonpara_logist_visit1.RData")
usethis::use_data(nonpara_logist_visit1, overwrite = TRUE)
