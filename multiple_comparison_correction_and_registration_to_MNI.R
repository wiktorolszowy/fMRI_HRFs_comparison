

###############################################################################################
####   Registering results to MNI space and doing FSL multiple testing. For AFNI, FSL and SPM.
####   Written by:    Wiktor Olszowy, University of Cambridge
####   Contact:       wo222@cam.ac.uk
####   Created:       September 2018 - November 2018
####   Adapted from:  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/multiple_comparison_correction_and_registration_to_MNI.R
###############################################################################################


#first run:           module load AFNI/AFNI_18.0.11/

library(reshape2)
library(tools)
library(R.matlab)
library(AnalyzeFMRI)
library(stringr)      #-for str_replace_all

path_manage        = readLines("path_manage.txt")
path_scratch       = readLines("path_scratch.txt")
studies_parameters = read.table(paste0(path_manage, "/studies_parameters.txt"), sep=";", header=T)
study              = paste0(studies_parameters$study[study_id])
abbr               = str_replace_all(paste0(studies_parameters$abbr [study_id]), pattern=" ", repl="")
task               = str_replace_all(paste0(studies_parameters$task [study_id]), pattern=" ", repl="")
packages           = c("AFNI", "FSL", "SPM")
exper_designs      = c(paste0("boxcar", c(12, 16, 20)), "event1", "event2")
subject            = paste0("sub-", abbr, strrep("0", 4-nchar(paste(subject_id))), subject_id)
path_data          = paste0(path_scratch, "/scans_", study)
path_output_top    = paste0(path_scratch, "/analysis_output_", study)

#-disabling scientific notation, otherwise an error for a volume of exactly 300000...
options(scipen = 999)


#-creating brain mask (intersected masks from FSL and SPM)
setwd(paste0(path_output_top, "/FSL/exper_design_boxcar12/HRF_gamma2_T/", subject))
system("fslmaths stats/zstat1 -abs zstat1_abs")
system(paste0("cp zstat1_abs.nii.gz ", path_output_top, "/together_masks/mask_FSL_", subject_id, ".nii.gz"))
setwd(paste0(path_output_top, "/SPM/exper_design_boxcar12/HRF_gamma2_T/", subject))
system(paste0("cp mask.nii ", path_output_top, "/together_masks/mask_SPM_", subject_id, ".nii"))
setwd(paste0(path_output_top, "/together_masks"))
system(paste0("rm mask_", subject_id, ".nii"))
system(paste0("fslmaths mask_FSL_", subject_id, " -bin mask_FSL_", subject_id))
system(paste0("fslmaths mask_FSL_", subject_id, " -add mask_SPM_", subject_id, " mask_", subject_id))
system(paste0("fslmaths mask_",     subject_id, " -thr 2 mask_", subject_id))
system(paste0("rm mask_FSL_", subject_id, ".nii.gz"))
system(paste0("rm mask_SPM_", subject_id, ".nii"))
system(paste0("gunzip mask_", subject_id, ".nii.gz"))

for (package in packages) {

   if (package=="AFNI") {
      HRF_models = c("gamma2", "gamma2_T", "gamma2_TD", "tent", "csplin")
   } else if (package=="FSL") {
      HRF_models = c("gamma2", "gamma2_T", "gamma", "gamma_T", "FIR")
   } else {
      HRF_models = c("gamma2", "gamma2_T", "gamma2_TD", "Fourier", "FIR", "gamma_x1", "gamma_x2", "gamma_x3")
   }

   for (exper_design in exper_designs) {

      for (HRF_model in HRF_models) {

         path_software = paste0(path_output_top, "/", package, "/exper_design_", exper_design, "/HRF_", HRF_model, "/", subject)

         setwd(path_software)

         system("rm -r standardized_stats")
         system("mkdir standardized_stats")

         cat(package, " ", exper_design, " ", HRF_model, "\n")

         if (length(list.files()) < 5 ) {
            cat("\n \n little output, probably the omnibus contrast from SPM/FAST did not return any voxels for BMMR!, at: \n", getwd(), "\n \n \n")
            next
         }

         if (HRF_model=="gamma" || HRF_model=="gamma2") {
            no_par = 1
         } else if (HRF_model=="gamma_T" || HRF_model=="gamma2_T") {
            no_par = 2
         } else if (HRF_model=="gamma2_TD") {
            no_par = 3
         } else if (HRF_model=="Fourier") {
            no_par = 11
         } else if (HRF_model=="FIR" || HRF_model=="tent" || HRF_model=="csplin") {
            no_par = 6
         } else if (HRF_model=="gamma_x1") {
            no_par = 1
         } else if (HRF_model=="gamma_x2") {
            no_par = 2
         } else if (HRF_model=="gamma_x3") {
            no_par = 3
         }
         
         if (package=="AFNI") {

            #-checking what are the degrees of freedom (df) of the t-statistic map
            #-https://afni.nimh.nih.gov/afni/community/board/read.php?1,67177,67179#msg-67179

            system(paste0("3dinfo stats.", subject, "_REML+orig['Full_Fstat'] > standardized_stats/3dinfo_output.txt"))
            info_output = readLines("standardized_stats/3dinfo_output.txt")
            line_20     = info_output[20]
            where_eq    = unlist(gregexpr(pattern="=", line_20))[2]
            dfs         = substr(line_20, where_eq+2, nchar(line_20))
            where_sp    = unlist(gregexpr(pattern=" ", dfs))
            df1         = no_par       #as.numeric(substr(dfs, 1, where_sp-1))
            df2         = as.numeric(substr(dfs, where_sp+1, nchar(dfs)))
            if (!is.numeric(df1) || !(df1>0) || !(df1<20) || !is.numeric(df2) || !(df2>10) || !(df2<1000)) {
               cat(paste0("FATAL ERROR related to df at ", getwd()))
            }

            #-transforming the F-statistic map to a z-statistic map
            #-https://afni.nimh.nih.gov/afni/community/board/read.php?1,157604,157604#msg-157604
            system(paste0("3dcalc -a stats.", subject, "_REML+orig'['Full_Fstat']' -expr 'cdf2stat(stat2cdf(a, 4, ", df1, ", ", df2, ", 0), 5, 0, 0, 0)' -prefix standardized_stats/zstat_with_extremes.nii"))

            #-there are many zeroes in the F-statistic map (primarily outside of the brain), which are transformed to extremely low numbers in zstat (e.g. -9.9900e+37)
            #-also, a lot of very high values in the F-statistic map, which are transformed to extremely high numbers in zstat
            system("fslmaths standardized_stats/zstat_with_extremes             -max -10000 standardized_stats/zstat_with_extremes_on_one_side")
            system("fslmaths standardized_stats/zstat_with_extremes_on_one_side -min  10000 standardized_stats/zstat")
            system("rm standardized_stats/zstat_with_extremes_on_one_side.nii.gz")
            
            #-transforming AFNI's GLM residuals to nifti
            #-ifs used to catch numerical problems with AFNI for 4 subjects in the "NKI RS TR=1.4s" dataset: ids 2&5 (no 'stats_REML'), 25&27 (no 'stats'):
            if (file.exists(paste0(    "whitened_errts.", subject, "_REML+orig.BRIK"))==T) {
               system(paste0('3dcalc -a whitened_errts.', subject, '_REML+orig -expr "a" -prefix standardized_stats/res4d.nii'));
            } else if (file.exists(paste0("errts.",       subject,      "+orig.BRIK"))==T) {
               system(paste0('3dcalc -a    errts.',       subject,      '+orig -expr "a" -prefix standardized_stats/res4d.nii'));
               cat(paste0("possible AFNI ERROR/problem (no 'whitened_errts' file) at ", getwd(), "\n"))
            } else {
               cat(paste0("possible AFNI ERROR/problem (no 'errts' file) at ", getwd(), "\n"))
               next
            }
            system("gzip standardized_stats/res4d.nii")
            
            #-moving 'coef's
            for (par_id in 1:no_par) {
               system(paste0("3dAFNItoNIFTI -prefix ", path_software, "/standardized_stats/coef_", par_id, " ", path_software, "/stats.", subject, "_REML+orig['activation_stimulus#", par_id-1, "_Coef']"), ignore.stderr=T)
            }
         
         } else if (package=="FSL") {

            #-for models with only 1 EV, FSL does not produce 'zfstat1'
            if (HRF_model=="gamma" || HRF_model=="gamma2") {
               system("cp stats/zstat1.nii.gz  standardized_stats/zstat.nii.gz")
            } else {
               system("cp stats/zfstat1.nii.gz standardized_stats/zstat.nii.gz")
            }

            #-reading the degrees of freedom
            df1 = no_par
            df2 = scan("stats/dof", quiet=T)

            #-compressing and copying FSL's GLM residuals
            if (file.exists("stats/res4d.nii")==T) {
               system("gzip stats/res4d.nii")
            }
            system("cp stats/res4d.nii.gz standardized_stats/res4d.nii.gz")

            #-moving 'coef's
            for (par_id in 1:no_par) {
               system(paste0("cp stats/pe", par_id, ".nii.gz standardized_stats/coef_", par_id, ".nii.gz"))
            }

         } else if (package=="SPM") {

            #-checking what are the degrees of freedom (df)
            df1 = no_par
            #-https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind06&L=SPM&P=R186259&1=SPM&9=A&I=-3&J=on&X=03B60AE4F53D3D3C56&Y=WO222%40CAM.AC.UK&d=No+Match%3BMatch%3BMatches&z=4
            #-SPM.xX.erdf
            SPM = readMat("SPM.mat")
            df2 = as.numeric(SPM$SPM[,,1]$xX[,,1]$erdf)

            #-transforming the F-statistic map to a z-statistic map
            #-https://afni.nimh.nih.gov/afni/community/board/read.php?1,157604,157604#msg-157604
            system(paste0("3dcalc -a spmF_0001.nii -expr 'cdf2stat(stat2cdf(a, 4, ", df1, ", ", df2, ", 0), 5, 0, 0, 0)' -prefix standardized_stats/zstat_with_extremes.nii"))

            #-there are many zeroes in the F-statistic map (primarily outside of the brain), which are transformed to extremely low numbers in zstat (e.g. -9.9900e+37)
            #-also, a lot of very high values in the F-statistic map, which are transformed to extremely high numbers in zstat
            system("fslmaths standardized_stats/zstat_with_extremes             -max -10000 standardized_stats/zstat_with_extremes_on_one_side")
            system("fslmaths standardized_stats/zstat_with_extremes_on_one_side -min  10000 standardized_stats/zstat")
            system("rm standardized_stats/zstat_with_extremes_on_one_side.nii.gz")

            #-merging SPM's GLM residuals to one file
            res_all = sort(list.files()[which(substr(list.files(), 1, 4)=="Res_")])
            res_all = paste(res_all, collapse=" ")
            system(paste0("fslmerge -t standardized_stats/res4d ", res_all))

            #-moving 'coef's
            for (par_id in 1:no_par) {
               if (par_id < 10) {
                  system(paste0("cp beta_000", par_id, ".nii standardized_stats/coef_", par_id, ".nii"))
               } else {
                  system(paste0("cp beta_00",  par_id, ".nii standardized_stats/coef_", par_id, ".nii"))
               }
            }

         }

         path = paste0(path_software, "/standardized_stats")
         setwd(path)

         #-applying brain mask to 'zstat' and 'res4d'
         system(paste0("cp ", path_output_top, "/together_masks/mask_", subject_id, ".nii mask.nii"))
         system(paste0("fslmaths ", path, "/zstat -mas mask ", path, "/zstat_masked"))
         system(paste0("fslmaths ", path, "/res4d -mas mask ", path, "/res4d_masked"))

         #-estimating smoothness, needed to run the 'cluster' function in FSL (multiple comparison correction)
         #-https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1710&L=FSL&D=0&X=CA0D5B74DD2B36BED5&Y=wo222%40cam.ac.uk&P=87052
         system(paste0("smoothest -r res4d_masked -d ", df2, " -m mask > smoothness"))
         par = read.table("smoothness", header=F, sep=" ", nrows=3)
         system(paste0("cluster --in=zstat_masked --oindex=cluster_F_index --thresh=3.09 --pthresh=0.05 --dlh=", par$V2[1], " --volume=", round(par$V2[2]), " > cluster_F.txt"))

         #-copying FSL MNI registration template and the corresponding transformation
         system(paste0("cp ", path_output_top, "/FSL/preproc_feats/", subject, "_preproc.feat/reg/standard.nii.gz standard.nii.gz"))
         if (study=="BMMR_checkerboard") {
            system(paste0("cp ", path_manage, "/BMMR_checkerboard_MNI_registrations/example_func2standard_", subject_id, ".mat example_func2standard.mat"))
         } else {
            system(paste0("cp ", path_output_top, "/FSL/preproc_feats/", subject, "_preproc.feat/reg/example_func2standard.mat example_func2standard.mat"))
         }

         #-changing all names of clusters to 1, to avoid later problems with interpolation
         system("fslmaths 'cluster_F_index' -thr 0.5 -bin 'cluster_binary_F'")

         #-registering 'cluster_binary_F' and 'mask' to MNI space
         system("flirt -ref standard -in cluster_binary_F -applyxfm -init example_func2standard.mat -out cluster_binary_F_MNI")
         system("flirt -ref standard -in mask             -applyxfm -init example_func2standard.mat -out mask_MNI")

         #-binarizing (needed due to MNI registration, which used interpolation)
         system("fslmaths 'cluster_binary_F_MNI' -thr 0.5 -bin 'cluster_binary_F_MNI'")
         system("fslmaths 'mask_MNI'             -thr 0.5 -bin 'mask_MNI'")

         #-applying MNI template mask (additional masking)
         system("fslmaths cluster_binary_F_MNI -mas standard cluster_binary_F_MNI")
         system("fslmaths mask_MNI             -mas standard mask_MNI")

         #-registering 'coef's to MNI space and masking
         for (par_id in 1:no_par) {
            system(paste0("flirt -ref standard -in coef_", par_id, " -applyxfm -init example_func2standard.mat -out coef_MNI_", par_id))
            system(paste0("fslmaths coef_MNI_", par_id, " -mas mask_MNI.nii coef_MNI_masked_", par_id))
            system(paste0("gunzip coef_MNI_masked_", par_id, ".nii"))
         }
         
         #-unzipping, because 'f.read.nifti.volume' does not work with '.gz'
         system("gunzip mask_MNI.nii.gz             -f")
         system("gunzip zstat_masked.nii.gz         -f")
         system("gunzip cluster_binary_F_MNI.nii.gz -f")
         mask                 = f.read.nifti.volume("mask.nii")
         mask_MNI             = f.read.nifti.volume("mask_MNI.nii")
         zstat_masked         = f.read.nifti.volume("zstat_masked.nii")
         cluster_binary_F_MNI = f.read.nifti.volume("cluster_binary_F_MNI.nii")

         #-0.5 chosen to distinguish clusters/brain voxels (>=1) from non-clusters/non-brain voxels (=0)
         mask_MNI             = mask_MNI > 0.5
         mask_MNI             = mask_MNI[,,,1] * 1.0
         cluster_binary_F_MNI = cluster_binary_F_MNI > 0.5
         cluster_binary_F_MNI = cluster_binary_F_MNI[,,,1] * 1.0

         #-'cluster_binary_F_MNI.mat' is the ultimate output, 3-dim array showing where the significant clusters are
         writeMat("cluster_binary_F_MNI.mat", cluster_binary_F_MNI = cluster_binary_F_MNI)

         #-saving the numbers of brain mask voxels and of significant voxels
         write(sum(mask_MNI             > 0.5), "no_of_MNI_mask_voxels")
         write(sum(cluster_binary_F_MNI > 0.5), "no_of_MNI_sig_voxels")

         #-checking how many voxels in 'zstat_masked' above 3.1
         zstat_masked_thr = zstat_masked     > 3.1
         zstat_masked_thr = zstat_masked_thr * 1.0

         #-saving the proportion of voxels with 'zstat_masked' above 3.1
         write(sum(zstat_masked_thr > 0.5) / sum(mask > 0.5), "prop_above_3_1")

         #-saving the size of 'res4d_masked'
         write(file.size("res4d_masked.nii.gz"), "res4d_masked_size")

         #-calculating and saving a 3D smoothness estimate (geometric mean of FWHMmm in x-, y- and z-directions)
         FWHMmm           = readLines("smoothness")
         FWHMmm           = FWHMmm[5]
         white_spaces     = gregexpr(' ', FWHMmm)[[1]][1:3]
         FWHMmm_x         = as.numeric(substr(FWHMmm, white_spaces[1]+1, white_spaces[2]-1))
         FWHMmm_y         = as.numeric(substr(FWHMmm, white_spaces[2]+1, white_spaces[3]-1))
         FWHMmm_z         = as.numeric(substr(FWHMmm, white_spaces[3]+1, nchar(FWHMmm)))
         FWHMmm_geom_mean = (FWHMmm_x*FWHMmm_y*FWHMmm_z)^(1/3)
         write(FWHMmm_geom_mean, "smoothness_3D")

         ##########################################################################################
         #-two-sided t-test

         if (HRF_model=="gamma2" || HRF_model=="gamma2_T" || HRF_model=="gamma2_TD") {

            setwd(path_software)
            
            if (package=="AFNI") {

               #-transforming the t-statistic map to a z-statistic map
               #-https://afni.nimh.nih.gov/afni/community/board/read.php?1,156394,156428#msg-156428
               system(paste0("3dcalc -a stats.", subject, "_REML+orig'['activation_stimulus#0_Tstat']' -expr 'fitt_t2z(a,", df2, ")' -prefix standardized_stats/zstat_t.nii"))

            } else if (package=="FSL") {

               system("cp stats/zstat1.nii.gz standardized_stats/zstat_t.nii.gz")

            } else if (package=="SPM") {

               system("cp spmT_0002.nii standardized_stats/spmT_0002.nii")

               #-'spmT_0002' rather than 'spmT_0001' was produced by SPM, as this was the second test (after F-test, which is why there is 'spmF_0001.nii')
               system("fslmaths standardized_stats/spmT_0002.nii -mul 0 -add 1 standardized_stats/ones.nii")
               system(paste0("ttoz standardized_stats/ones standardized_stats/spmT_0002 ", df2, " -zout standardized_stats/zstat_t"))

            }
            
            #-move to within the 'standardized_stats' folder
            setwd(path)
            
            #-applying brain mask (intersected from FSL and SPM brain masks)
            system(paste0("fslmaths ", path, "/zstat_t -mas ", path, "/mask ", path, "/zstat_t_masked"))

            #-absolute values taken (needed for two-sided t-test)
            system(paste0("fslmaths ", path, "/zstat_t_masked -abs ",          path, "/zstat_t_masked_abs"))

            #-cluster inference; as now absolute values used, significance level down from 5% to 2.5%
            system(paste0("cluster --in=zstat_t_masked_abs --oindex=cluster_t_index --thresh=3.09 --pthresh=0.025 --dlh=", par$V2[1], " --volume=", round(par$V2[2]), " > cluster_t.txt"))

            #-changing all names of clusters to 1, to avoid later problems with interpolation
            system("fslmaths 'cluster_t_index' -thr 0.5 -bin 'cluster_binary_t'")

            #-registering results to MNI space
            system("flirt -ref standard -in cluster_binary_t -applyxfm -init example_func2standard.mat -out cluster_binary_t_MNI")

            #-binarizing (needed due to MNI registration, which used interpolation)
            system("fslmaths 'cluster_binary_t_MNI' -thr 0.5 -bin 'cluster_binary_t_MNI'")

            #-applying MNI template mask (additional masking)
            system("fslmaths cluster_binary_t_MNI -mas standard cluster_binary_t_MNI")

            #-unzipping, because 'f.read.nifti.volume' does not work with '.gz'
            system("gunzip cluster_binary_t_MNI.nii.gz -f")
            cluster_binary_t_MNI = f.read.nifti.volume("cluster_binary_t_MNI.nii")

            cluster_binary_t_MNI = cluster_binary_t_MNI > 0.5
            cluster_binary_t_MNI = cluster_binary_t_MNI[,,,1] * 1.0

            #-'cluster_binary_t_MNI.mat' is the ultimate output, 3-dim array showing where the significant clusters are
            writeMat("cluster_binary_t_MNI.mat", cluster_binary_t_MNI = cluster_binary_t_MNI)

         }
         ##########################################################################################

         #-to make space
         system("rm res4d.nii.gz")

      }

   }

}

setwd(path_manage)
