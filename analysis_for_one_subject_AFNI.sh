#!/bin/bash

###############################################################################################
####   AFNI analysis for 1 fMRI scan, for different combinations of options.
####   Written by:    Wiktor Olszowy, University of Cambridge
####   Contact:       wo222@cam.ac.uk
####   Created:       September 2018 - November 2018
####   Adapted from:  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/analysis_for_one_subject_AFNI.sh
###############################################################################################


path_manage=`cat path_manage.txt`
path_scratch=`cat path_scratch.txt`
cd $path_manage
declare -a studies=$(awk -F';' '{ if (NR!=1) print $1 }' studies_parameters.txt)
declare -a     TRs=$(awk -F';' '{ if (NR!=1) print $2 }' studies_parameters.txt)
declare -a      ns=$(awk -F';' '{ if (NR!=1) print $5 }' studies_parameters.txt)
declare -a   abbrs=$(awk -F';' '{ if (NR!=1) print $6 }' studies_parameters.txt)
declare -a   tasks=$(awk -F';' '{ if (NR!=1) print $7 }' studies_parameters.txt)
studies=($studies)
TRs=($TRs)
ns=($ns)
abbrs=($abbrs)
tasks=($tasks)
study=${studies[study_id-1]}
TR=${TRs[study_id-1]}
n=${ns[study_id-1]}
abbr=${abbrs[study_id-1]}
task=${tasks[study_id-1]}
no_zeros=`echo 4-${#subject_id} | bc`
if [ ${no_zeros} -ge 1 ]; then
   zeros="$(printf '0%.0s' $(seq 1 $no_zeros))"
else
   zeros=""
fi
printf "$1"'%.s' $(eval "echo {1.."$(($2))"}");
subject=sub-${abbr}$zeros${subject_id}
declare -a smoothing=5
declare -a exper_designs=("boxcar12" "boxcar16" "boxcar20" "event1" "event2")
declare -a HRF_models=("gamma2" "gamma2_T" "gamma2_TD" "tent" "csplin")
path_data=${path_scratch}/scans_${study}
path_output=${path_scratch}/analysis_output_${study}/AFNI

module load AFNI/AFNI_18.0.11/

#-defining slice timing correction
case $study in
   CamCAN_sensorimotor)
      tpattern=seqminus
      ;;
   BMMR_checkerboard)
      tpattern=alt+z
      ;;
   CRIC_checkerboard)
      tpattern=alt+z2
      ;;
esac

for ((exper_design_id=0; exper_design_id<${#exper_designs[@]}; exper_design_id++)) {

   exper_design=${exper_designs[exper_design_id]}

   for ((HRF_model_id=0; HRF_model_id<${#HRF_models[@]}; HRF_model_id++)) {

      #-defining the HRF approach
      HRF_model=${HRF_models[HRF_model_id]}
      #-stimulus duration in brackets needed for convolution
      if [ ${exper_design:0:6} == boxcar ]; then
         HRF_dur=${exper_design:6:2}
      else
         HRF_dur=0
      fi
      case $HRF_model in
         gamma2)
            HRF="SPMG1(${HRF_dur})"
            ;;
         gamma2_T)
            HRF="SPMG2(${HRF_dur})"
            ;;
         gamma2_TD)
            HRF="SPMG3(${HRF_dur})"
            ;;
         tent)
            HRF="TENT(3,18,6)"
            ;;
         csplin)
            HRF="CSPLIN(3,18,6)"
            ;;
      esac

      echo $study $subject $exper_design $HRF_model

      #-defining path to the stimulus times
      if [ ${study} == CamCAN_sensorimotor ] && [ ${exper_design} == event1 ]; then
         stim_times=${path_manage}/experimental_designs/CamCAN_sensorimotor/AFNI_$zeros${subject_id}.txt
      else
         stim_times=${path_manage}/experimental_designs/AFNI_${study}_${exper_design}.txt
      fi
      
      cd ${path_output}/exper_design_${exper_design}/HRF_${HRF_model}
      rm proc.${subject}
      rm output.proc.${subject}
      rm -rf ${subject}
      rm -rf ${subject}.results
      
      #-read https://afni.nimh.nih.gov/afni/community/board/read.php?1,156558,156558#msg-156558
      if [ "$study" == BMMR_checkerboard ]; then
         #-without motion correction
         #-with slice timing correction
         afni_proc.py                                                                \
            -subj_id $subject                                                        \
            -script proc.$subject -scr_overwrite                                     \
            -blocks tshift blur mask scale regress                                   \
            -tshift_opts_ts -tpattern $tpattern                                      \
            -dsets ${path_data}/${subject}_${task}_bold.nii                          \
            -regress_polort 2                                                        \
            -regress_bandpass 0.01 10                                                \
            -blur_size ${smoothing}                                                  \
            -regress_stim_types AM1                                                  \
            -regress_stim_times $stim_times                                          \
            -regress_stim_labels activation_stimulus                                 \
            -regress_basis ${HRF}                                                    \
            -regress_make_ideal_sum sum_ideal.1D                                     \
            -regress_run_clustsim no                                                 \
            -regress_est_blur_epits                                                  \
            -regress_est_blur_errts                                                  \
            -regress_reml_exec                                                       \
            -regress_opts_reml -Rwherr whitened_errts.${subject}_REML
      elif [ "$study" == NKI_release_3_checkerboard_1400 ] || [ "$study" == NKI_release_3_checkerboard_645 ] ; then
         #-with motion correction
         #-without slice timing correction
         afni_proc.py                                                                \
            -subj_id $subject                                                        \
            -script proc.$subject -scr_overwrite                                     \
            -dsets ${path_data}/${subject}_${task}_bold.nii                          \
            -volreg_align_to third                                                   \
            -regress_polort 2                                                        \
            -regress_bandpass 0.01 10                                                \
            -blur_size ${smoothing}                                                  \
            -regress_stim_types AM1                                                  \
            -regress_stim_times $stim_times                                          \
            -regress_stim_labels activation_stimulus                                 \
            -regress_basis ${HRF}                                                    \
            -regress_make_ideal_sum sum_ideal.1D                                     \
            -regress_run_clustsim no                                                 \
            -regress_est_blur_epits                                                  \
            -regress_est_blur_errts                                                  \
            -regress_reml_exec                                                       \
            -regress_opts_reml -Rwherr whitened_errts.${subject}_REML
      else
         #-with motion correction
         #-with slice timing correction
         afni_proc.py                                                                \
            -subj_id $subject                                                        \
            -script proc.$subject -scr_overwrite                                     \
            -tshift_opts_ts -tpattern $tpattern                                      \
            -dsets ${path_data}/${subject}_${task}_bold.nii                          \
            -volreg_align_to third                                                   \
            -regress_polort 2                                                        \
            -regress_bandpass 0.01 10                                                \
            -blur_size ${smoothing}                                                  \
            -regress_stim_types AM1                                                  \
            -regress_stim_times $stim_times                                          \
            -regress_stim_labels activation_stimulus                                 \
            -regress_basis ${HRF}                                                    \
            -regress_make_ideal_sum sum_ideal.1D                                     \
            -regress_run_clustsim no                                                 \
            -regress_est_blur_epits                                                  \
            -regress_est_blur_errts                                                  \
            -regress_reml_exec                                                       \
            -regress_opts_reml -Rwherr whitened_errts.${subject}_REML
      fi

      #-run the AFNI analysis
      tcsh -xef proc.${subject} |& tee  output.proc.${subject}
      
      mv ${subject}.results ${subject}
      mv proc.${subject} ${subject}/proc.${subject}
      mv output.proc.${subject} ${subject}/output.proc.${subject}

   }

}

cd ${path_manage}
