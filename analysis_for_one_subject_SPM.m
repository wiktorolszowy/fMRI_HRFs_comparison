

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   SPM analysis for 1 fMRI scan.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       September 2018 - November 2018
%%%%   Adapted from:  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/analysis_for_one_subject_SPM.m
%%%%                  https://github.com/wanderine/ParametricMultisubjectfMRI/blob/master/SPM/analyze_all_subjects_spm_1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage        = fgetl(fopen('path_manage.txt'));
path_scratch       = fgetl(fopen('path_scratch.txt'));
studies_parameters = readtable([path_manage '/studies_parameters.txt']);
smoothing          = 5;
exper_designs      = cellstr(['boxcar12'; 'boxcar16'; 'boxcar20'; 'event1  '; 'event2  ']);
HRF_models         = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 ']);
stim_onsets_CamCAN = textscan(fopen([path_manage '/experimental_designs/CamCAN_sensorimotor/SPM_' repmat('0', 1, 4-length(num2str(subject_id))) num2str(subject_id) '.txt']), '%f');
stim_onsets_event1 = textscan(fopen([path_manage '/experimental_designs/SPM_event1.txt']), '%f');
stim_onsets_event2 = textscan(fopen([path_manage '/experimental_designs/SPM_event2.txt']), '%f');
study              = studies_parameters.study{study_id};
abbr               = studies_parameters.abbr{study_id};
task               = studies_parameters.task{study_id};
subject            = ['sub-' abbr repmat('0', 1, 4-length(num2str(subject_id))) num2str(subject_id)];
bold_file          = [subject '_' task '_bold'];
TR                 = studies_parameters.TR(study_id);                      %-the repetition time
npts               = studies_parameters.npts(study_id);                    %-the number of time points
max_index_event1   = find((stim_onsets_event1{1}+4)<npts*TR, 1, 'last');
max_index_event2   = find((stim_onsets_event2{1}+4)<npts*TR, 1, 'last');
path_data          = [path_scratch '/scans_' study];
path_preproc_top   = [path_scratch '/analysis_output_' study '/SPM/preproc'];
path_preproc       = [path_preproc_top '/' subject];
filename           = [path_preproc '/' bold_file '.nii'];

addpath('/applications/spm/spm12_7219');
addpath(genpath([path_manage '/matlab_extra_functions']));

%-initialise SPM defaults
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clear jobs;

%-default in SPM is 128, but in FSL it is 100
spm_get_defaults('stats.fmri.hpf', 100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(['rm -rf ' path_preproc_top '/' subject]);
system(['mkdir  ' path_preproc_top '/' subject]);

%-copy the scan
system(['cp ' path_data '/' bold_file '.nii '           path_preproc]);

cd(path_preproc)

%-REALIGN (motion correction)
if strcmp(study, 'BMMR_checkerboard')
   %-for 'BMMR_checkerboard', motion correction was applied before the data was distributed
   system(['mv ' bold_file '.nii r' bold_file '.nii']);
else
   clear jobs;
   jobs{1}.spatial{1}.realign{1}.estwrite.data{1} = cellstr(filename);
   spm_jobman('run', jobs);
end

%-slice timing correction
%-page 248 of https://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf
clear jobs;
if strcmp(study, 'CamCAN_sensorimotor')
   nslices     = 32;
   slice_order = nslices:-1:1;
elseif strcmp(study, 'NKI_release_3_checkerboard_1400')
   nslices     = 64;
   slice_order = round(nslices/2) * ones(nslices,1);
elseif strcmp(study, 'NKI_release_3_checkerboard_645')
   nslices     = 40;
   slice_order = round(nslices/2) * ones(nslices,1);
elseif strcmp(study, 'BMMR_checkerboard')
   nslices     = 45;
   slice_order = [1:2:nslices 2:2:nslices];
elseif strcmp(study, 'CRIC_checkerboard')
   nslices     = 32;
   slice_order = [2:2:nslices 1:2:nslices];
end
jobs{1}.temporal{1}.st.scans{1} = spm_select('expand', cellstr(spm_file(filename, 'prefix', 'r')));
jobs{1}.temporal{1}.st.nslices  = nslices;
jobs{1}.temporal{1}.st.tr       = TR;
jobs{1}.temporal{1}.st.ta       = TR - (TR/nslices);
jobs{1}.temporal{1}.st.so       = slice_order;
jobs{1}.temporal{1}.st.refslice = round(nslices/2);

%-SMOOTHING
jobs{2}.spatial{1}.smooth.data  = cellstr(spm_file(filename, 'prefix', 'ar'));
jobs{2}.spatial{1}.smooth.fwhm  = [smoothing smoothing smoothing];

spm_jobman('run', jobs);

%-registration to MNI space is later conducted using transformations from FSL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-STATISTICAL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for exper_design_id = 1:length(exper_designs)

   exper_design     = exper_designs{exper_design_id};

   for HRF_model_id = 1:5

      HRF_model     = HRF_models{HRF_model_id};

      %-to load the appropriate 'spm_get_bf.m'; that is where the combined gamma functions (modelling main response, undershoot and initial dip) are implemented
      if strcmp(HRF_model, 'gamma_x1') || strcmp(HRF_model, 'gamma_x2') || strcmp(HRF_model, 'gamma_x3')
         addpath(genpath([path_manage '/matlab_extra_functions']));
      else
         addpath('/applications/spm/spm12_7219');
      end

      path_output   = [path_scratch '/analysis_output_' study '/SPM/exper_design_' exper_design '/HRF_' HRF_model];
      cd(path_output);
      system(['rm -rf ' subject]);
      system(['rm -rf ' subject '_old']);
      system(['mkdir ' subject]);

      if (study_id==4) && (subject_id==13) && (exper_design_id==1) && (HRF_model_id==1)
         %-some weird error
         continue
      end
      
      disp([study, ' ', exper_design, ' ', HRF_model, ' ', subject]);

      clear jobs;
      
      %-MODEL SPECIFICATION AND ESTIMATION
      filename                                = [path_output '/' subject];
      jobs{1}.stats{1}.fmri_spec.dir          = cellstr(filename);
      jobs{1}.stats{1}.fmri_spec.timing.units = 'secs';
      jobs{1}.stats{1}.fmri_spec.timing.RT    = TR;
      
      scans = {};
      for t = 1:npts
         scans{t} = [path_preproc '/sar' bold_file '.nii,' num2str(t)];
      end
      jobs{1}.stats{1}.fmri_spec.sess.scans        = transpose(scans);
      jobs{1}.stats{1}.fmri_spec.sess.cond(1).name = 'task1';

      %-specifying the experimental design
      if strcmp(exper_design(1:6), 'boxcar')
         len = str2num(exper_design(7:8));
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = len:(2*len):(npts*TR-len);
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = len;
      elseif strcmp(study, 'CamCAN_sensorimotor') && strcmp(exper_design, 'event1')
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = stim_onsets_CamCAN{1};
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = 0.1;
      elseif strcmp(exper_design, 'event1')
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = stim_onsets_event1{1}(1:max_index_event1);
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = 0.1;
      elseif strcmp(exper_design, 'event2')
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = stim_onsets_event2{1}(1:max_index_event2);
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = 0.1;
      end

      %-SPM's default is 0.8, but because some masks for 0.8 were too small, I chose a more liberal threshold
      jobs{1}.stats{1}.fmri_spec.mthresh                  = 0.2;
      
      %-including motion regressors in the GLM
      if ~strcmp(study, 'BMMR_checkerboard')
         jobs{1}.stats{1}.fmri_spec.sess.multi_reg        = {[path_preproc '/rp_' bold_file '.txt']};
      end

      if strcmp(HRF_model, 'gamma2')
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs      = [0 0];
         no_par                                           = 1;
      elseif strcmp(HRF_model, 'gamma2_T')
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs      = [1 0];
         no_par                                           = 2;
      elseif strcmp(HRF_model, 'gamma2_TD')
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs      = [1 1];
         no_par                                           = 3;
      elseif strcmp(HRF_model, 'Fourier')
         jobs{1}.stats{1}.fmri_spec.bases.fourier.length  = 12*1.97;
         jobs{1}.stats{1}.fmri_spec.bases.fourier.order   = 5;
         no_par                                           = 11;
      elseif strcmp(HRF_model, 'FIR')
         jobs{1}.stats{1}.fmri_spec.bases.fir.length      = 6*3;
         jobs{1}.stats{1}.fmri_spec.bases.fir.order       = 6;
         no_par                                           = 6;
      elseif strcmp(HRF_model, 'gamma_x1')
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs      = [0 0];
         no_par                                           = 1;
      elseif strcmp(HRF_model, 'gamma_x2')
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs      = [1 0];
         no_par                                           = 2;
      elseif strcmp(HRF_model, 'gamma_x3')
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs      = [1 1];
         no_par                                           = 3;
      end
      
      %-SPM's alternative pre-whitening: FAST
      %-following study "Accurate autocorrelation modeling substantially improves fMRI reliability":
      %-https://www.biorxiv.org/content/early/2018/09/02/323154
      jobs{1}.stats{1}.fmri_spec.cvi              = 'FAST';

      filename_mat                                = [filename '/SPM.mat'];
      jobs{2}.stats{1}.fmri_est.spmmat            = cellstr(filename_mat);

      %-F-test
      jobs{3}.stats{1}.con.spmmat                 = cellstr(filename_mat);
      jobs{3}.stats{1}.con.consess{1}.fcon        = struct('name', 'task1 > rest', 'convec', eye(no_par), 'sessrep', 'none');

      %-additionally t-test
      if strcmp(HRF_model, 'gamma2') || strcmp(HRF_model, 'gamma2_T') || strcmp(HRF_model, 'gamma2_TD')
         jobs{4}.stats{1}.con.spmmat              = cellstr(filename_mat);
         jobs{4}.stats{1}.con.consess{1}.tcon     = struct('name', 'task1 > rest', 'convec', 1,           'sessrep', 'none');
      end

      spm_jobman('run', jobs);

      %-saving/not removing GLM residuals
      cd([path_output '/' subject]);
      load('SPM.mat');
      VRes = spm_write_residuals(SPM, NaN);
      
   end

end

cd(path_manage)
