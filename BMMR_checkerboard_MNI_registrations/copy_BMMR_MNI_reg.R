
for (subject_id in 1:21) {
   if (subject_id < 10) {
      system(paste0("cp /home/wo222/scratch/fMRI_method_validations/FINISHED_modelling_autocorr/analysis_output_BMMR_checkerboard/AFNI/smoothing_8/exper_design_boxcar10/HRF_gamma2_T/sub-BMM000", subject_id, "/standardized_stats/example_func2standard.mat /home/wo222/fMRI_method_validations/modelling_HRF_comparison/analysis/BMMR_checkerboard_MNI_registrations/example_func2standard_", subject_id, ".mat"))
   } else {
      system(paste0("cp /home/wo222/scratch/fMRI_method_validations/FINISHED_modelling_autocorr/analysis_output_BMMR_checkerboard/AFNI/smoothing_8/exper_design_boxcar10/HRF_gamma2_T/sub-BMM00",  subject_id, "/standardized_stats/example_func2standard.mat /home/wo222/fMRI_method_validations/modelling_HRF_comparison/analysis/BMMR_checkerboard_MNI_registrations/example_func2standard_", subject_id, ".mat"))
   }
}
