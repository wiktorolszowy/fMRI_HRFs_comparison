
#-Wiktor Olszowy

library(parallel)

path_func    = "/home/wo222/scratch/fMRI_data/CamCAN/cc700/mri/pipeline/release004/BIDSsep/func_smt"
path_rest    = "/home/wo222/scratch/fMRI_data/CamCAN/cc700/mri/pipeline/release004/BIDSsep/func_rest"
path_anat    = "/home/wo222/scratch/fMRI_data/CamCAN/cc700/mri/pipeline/release004/BIDSsep/anat"
path_smt_to  = "/home/wo222/scratch/fMRI_method_validations/modelling_HRF/scans_CamCAN_sensorimotor"
path_rest_to = "/home/wo222/scratch/fMRI_method_validations/modelling_HRF/scans_CamCAN_RS"

setwd(path_func)
subjects = list.files()
subjects = subjects[2:length(subjects)]

for (subject in subjects) {
   setwd(paste0(path_func, "/", subject, "/func"))
   if (length(which(list.files()==paste0(subject, "_task-SMT_bold.json")))!=1) {
      cat(subject, " func json \n")
   }
   if (length(which(list.files()==paste0(subject, "_task-SMT_bold.nii.gz")))!=1) {
      cat(subject, " func nii.gz \n")
   }
   if (length(which(list.files()==paste0(subject, "_task-SMT_events.tsv")))!=1) {
      cat(subject, " func tsv \n")
   }
   setwd(paste0(path_anat, "/", subject, "/anat"))
   if (length(which(list.files()==paste0(subject, "_T1w.nii.gz")))!=1) {
      cat(subject, " anat \n")
   }
   setwd(paste0(path_rest, "/", subject, "/func"))
   if (length(which(list.files()==paste0(subject, "_task-Rest_bold.nii.gz")))!=1) {
      cat(subject, " rest \n")
   }
   if (length(which(list.files()==paste0(subject, "_task-Rest_bold.json")))!=1) {
      cat(subject, " rest json \n")
   }
}
#-remove subjects where at least one file is missing
subjects_out = c("sub-CC110045", "sub-CC120065", "sub-CC221755", "sub-CC410222", "sub-CC610146")
subjects = subjects[-which(subjects %in% subjects_out)]

subjects = sort(subjects)
#-remove subjects where one of the scans has unusual size
subjects = subjects[-c(16, 93, 194, 226, 300, 308, 339, 412, 418, 422, 455, 461, 471, 473, 481, 491, 504, 571, 577, 588, 595)]

res = mclapply(1:length(subjects), function(subject_id) {
   subject = subjects[subject_id]
   sub_name_new = paste0("sub-CAS", strrep("0", 4-nchar(paste(subject_id))), subject_id)
   system(paste0("cp ", path_func, "/", subject, "/func/", subject, "_task-SMT_bold.json ",    path_smt_to,  "/", sub_name_new, "_sensorimotor_bold.json"))
   system(paste0("cp ", path_func, "/", subject, "/func/", subject, "_task-SMT_bold.nii.gz ",  path_smt_to,  "/", sub_name_new, "_sensorimotor_bold.nii.gz"))
   system(paste0("cp ", path_func, "/", subject, "/func/", subject, "_task-SMT_events.tsv ",   path_smt_to,  "/", sub_name_new, "_sensorimotor_events.tsv"))
   system(paste0("cp ", path_anat, "/", subject, "/anat/", subject, "_T1w.nii.gz ",            path_smt_to,  "/", sub_name_new, "_T1w.nii.gz"))
   system(paste0("cp ", path_rest, "/", subject, "/func/", subject, "_task-Rest_bold.json ",   path_rest_to, "/", sub_name_new, "_rest_bold.json"))
   system(paste0("cp ", path_rest, "/", subject, "/func/", subject, "_task-Rest_bold.nii.gz ", path_rest_to, "/", sub_name_new, "_rest_bold.nii.gz"))
   system(paste0("cp ", path_anat, "/", subject, "/anat/", subject, "_T1w.nii.gz ",            path_rest_to, "/", sub_name_new, "_T1w.nii.gz"))
}, mc.cores=24)

for (path in c(path_smt_to, path_rest_to)) {
   setwd(path)
   files = list.files()
   compressed = files[which(substr(files, nchar(files)-5, nchar(files))=="nii.gz")]
   res = mclapply(compressed, function(file) {
      system(paste0("gunzip ", file))
   }, mc.cores=24)
}

setwd(path_rest_to)
files = list.files()
for (file in files) {
   system(paste0("mv ", file, " sub-CAR", substr(file, 8, nchar(file))))
}

#-extract brains
setwd(path_smt_to)
files = list.files()
scans_T1 = files[which(substr(files, nchar(files)-6, nchar(files))=="T1w.nii")]

res = mclapply(scans_T1, function(scan_T1) {
   system(paste0("bet ", scan_T1, " ", substr(scan_T1, 1, nchar(scan_T1)-4), "_brain.nii -f 0.25"))
}, mc.cores=24)

files = list.files()
compressed = files[which(substr(files, nchar(files)-5, nchar(files))=="nii.gz")]
res = mclapply(compressed, function(file) {
   system(paste0("gunzip ", file))
}, mc.cores=24)


library(parallel)

setwd("/home/wo222/scratch/fMRI_method_validations/modelling_HRF_comparison/scans_CamCAN_sensorimotor")
scans = list.files()
scans = scans[which(substr(scans, nchar(scans)-8, nchar(scans))=="_bold.nii")]

res = mclapply(1:621, function(scan_id) {
   system(paste0("3drefit -TR 1.97 -view orig -space ORIG ", scans[scan_id]))
}, mc.cores=24)
