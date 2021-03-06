
#-Wiktor Olszowy

studies_parameters = read.table("studies_parameters.txt", sep=";", header=T)
subject_nos        = studies_parameters$n
lims               = seq(1, sum(subject_nos)*9+1, by=sum(subject_nos))
part               = 1

system("mkdir parallel_commands")

#-AFNI
id                 = 0
for (i in 1:length(subject_nos)) {
   for (j in 1:subject_nos[i]) {
      id = id + 1
      cat("export study_id=", i, "; export subject_id=", j, "; bash analysis_for_one_subject_AFNI.sh; \n", file=paste0("parallel_commands/command_", part, "_", id, ".sh"), sep="", append=F)
   }
}

#-FSL
id     = 0
part   = part+1
for (i in 1:length(subject_nos)) {
   for (j in 1:subject_nos[i]) {
      id = id + 1
      cat("R -e  'study_id=", i, "; subject_id=", j, "; source(\"analysis_for_one_subject_FSL.R\")' \n", file=paste0("parallel_commands/command_", part,   "_", id, ".sh"), sep="", append=F)
   }
}

#-SPM
id   = 0
part = part+1
for (i in 1:length(subject_nos)) {
   for (j in 1:subject_nos[i]) {
      id = id + 1
      cat("matlab -r -nodesktop \"study_id=", i, "; subject_id=", j, "; run('analysis_for_one_subject_SPM.m'); exit\" \n", file=paste0("parallel_commands/command_", part,   "_", id, ".sh"), sep="", append=F)
   }
}

#-perform multiple comparison correction and register results to MNI space
id     = 0
part   = part+1
for (i in 1:length(subject_nos)) {
   for (j in 1:subject_nos[i]) {
      id = id + 1
      cat("R -e  'study_id=", i, "; subject_id=", j, "; source(\"multiple_comparison_correction_and_registration_to_MNI.R\")' \n", file=paste0("parallel_commands/command_", part,   "_", id, ".sh"), sep="", append=F)
   }
}
