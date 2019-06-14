## logist model

load("C:/Users/unkyunglee/Dropbox/Project_Posdoc_Tanya/HDproject_local/Rcode(Revision_Rpackage)/code_for_R_package/HDChangePoint/Data/nlme_logist_visit1.RData")

nlme_logist_visit1<-results
save(nlme_logist_visit1, file="Data/nlme_logist_visit1.RData")
usethis::use_data(nlme_logist_visit1, overwrite = TRUE)


## arctangent model
load("C:/Users/unkyunglee/Dropbox/Project_Posdoc_Tanya/HDproject_local/Rcode(Revision_Rpackage)/code_for_R_package/HDChangePoint/Data/nlme_arctan_visit1.RData")

nlme_arctan_visit1<-results
save(nlme_arctan_visit1, file="Data/nlme_arctan_visit1.RData")
usethis::use_data(nlme_arctan_visit1, overwrite = TRUE)


## misspecified arctangent model
load("C:/Users/unkyung/Dropbox/Project_Posdoc_Tanya/HDproject_local/Rcode(Revision_Rpackage)/code_for_R_package/HDChangePoint/Data/arctan_missp_visit1.RData")

arctan_missp_visit1<-results
save(arctan_missp_visit1, file="Data/arctan_missp_visit1.RData")
usethis::use_data(arctan_missp_visit1, overwrite = TRUE)



## misspecified logist model

load("C:/Users/unkyung/Dropbox/Project_Posdoc_Tanya/HDproject_local/Rcode(Revision_Rpackage)/code_for_R_package/HDChangePoint/Data/missp_logist_visit1.RData")

missp_logist_visit1<-results
save(missp_logist_visit1, file="Data/missp_logist_visit1.RData")
usethis::use_data(missp_logist_visit1, overwrite = TRUE)

