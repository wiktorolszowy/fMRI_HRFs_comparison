

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Calculating power spectra on GLM residuals from AFNI/FSL/SPM.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       September 2018 - November 2018
%%%%   Adapted from:  https://github.com/wanderine/ParametricSinglesubjectfMRI/blob/master/FSL/fsl_powerspectra.m
%%%%                  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/make_power_spectra.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage         = fgetl(fopen('path_manage.txt'));
path_scratch        = fgetl(fopen('path_scratch.txt'));
path_output         = [path_scratch '/analysis_output_'];
studies_parameters  = readtable([path_manage '/studies_parameters.txt']);
studies             = studies_parameters.study;
packages            = cellstr(['AFNI'; 'FSL '; 'SPM ']);
exper_designs       = cellstr(['boxcar12'; 'boxcar16'; 'boxcar20'; 'event1  '; 'event2  ']);
p                   = parpool(12);
range_studies       = 1:length(studies);
range_packages      = 1:length(packages);
range_exper_designs = 1:length(exper_designs);
fft_n               = 512; %FFT will pad the voxel-wise time series to that length with trailing zeros (if no. of time points lower) or truncate to that length (if no. of time points higher)


cd(path_manage);
addpath(genpath([path_manage '/matlab_extra_functions']));
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');


for study_id   = range_studies

   study       = studies_parameters.study{study_id};
   abbr        = studies_parameters.abbr{study_id};
   task        = studies_parameters.task{study_id};
   no_subjects = studies_parameters.n(study_id);

   for package_id = range_packages

      package     = packages{package_id};

      if strcmp(package, 'AFNI')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'tent     '; 'csplin   ']);
      elseif strcmp(package, 'FSL')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma    '; 'gamma_T  '; 'FIR      ']);
      elseif strcmp(package, 'SPM')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 ']);
      end

      for exper_design_id = range_exper_designs

         exper_design     = exper_designs{exper_design_id};

         for HRF_model_id = 1:length(HRF_models)

            HRF_model     = HRF_models{HRF_model_id};

            disp([study ' ' package ' ' exper_design ' ' HRF_model]);

            path_subjects     = [path_output study '/' package '/exper_design_' exper_design '/HRF_' HRF_model];
            cd(path_subjects);

            %system('mv power_spectra.mat power_spectra_old.mat');
            if exist('power_spectra.mat', 'file') == 2
               system('rm power_spectra.mat');
            end

            parfor subject_id = 1:no_subjects

               subject        = ['sub-' abbr repmat('0', 1, 4-length(num2str(subject_id))) num2str(subject_id)];

               if exist([path_subjects '/' subject '/standardized_stats/cluster_binary_F_MNI.nii'], 'file') ~= 2
                  %-subjects without full output/without power spectra (some rare BMMR runs, only SPM/FAST)
                  continue
               end

               cd([path_subjects '/' subject '/standardized_stats']);

               %-removing output of previous runs
               if exist('power_spectra_one_subject.mat', 'file') == 2
                  system('rm power_spectra_one_subject.mat');
               end

               %-uncompressing, because niftiread has some rare/irreproducible numerical problems when reading nii.gz
               if exist('res4d_masked.nii.gz',  'file') == 2
                  if exist('res4d_masked.nii',  'file') == 2
                     system('mv res4d_masked.nii res4d_masked_old.nii');
                  end
                  system('gunzip res4d_masked.nii.gz');
               end
               if exist('zstat_masked.nii.gz', 'file') == 2
                  if exist('zstat_masked.nii', 'file') == 2
                     system('mv zstat_masked.nii zstat_masked_old.nii');
                  end
                  system('gunzip zstat_masked.nii.gz');
               end

               res4d        = niftiread('res4d_masked.nii');
               zstat_masked = niftiread('zstat_masked.nii');

               system('gzip res4d_masked.nii');
               system('gzip zstat_masked.nii');

               dims                      = size(res4d);
               power_spectra_one_subject = zeros(fft_n, 1);

               for i1 = 1:dims(1)

                  for i2 = 1:dims(2)

                     for i3 = 1:dims(3)

                        if zstat_masked(i1, i2, i3) ~= 0

                           ts = squeeze(res4d(i1, i2, i3, :));

                           if (std(ts) ~= 0)

                              %-make signal variance equal to 1
                              ts  = ts/(std(ts) + eps);

                              %-compute the discrete Fourier transform (DFT)
                              DFT = fft(ts, fft_n);

                              power_spectra_one_subject = power_spectra_one_subject + ((abs(DFT)).^2)/min(dims(4), fft_n);

                           end

                        end

                     end

                  end

               end

               %-average power spectra over all brain voxels
               power_spectra_one_subject = power_spectra_one_subject / sum(zstat_masked(:) ~= 0);

               %-trick to save within parfor
               parsave('power_spectra_one_subject.mat', power_spectra_one_subject);

            end

            power_spectra  = zeros(fft_n, 1);

            %-number of subjects without full output/without power spectra (some rare BMMR runs, only SPM/FAST)
            no_wo_output   = 0;

            for subject_id = 1:no_subjects

               subject     = ['sub-' abbr repmat('0', 1, 4-length(num2str(subject_id))) num2str(subject_id)];

               if exist([path_subjects '/' subject '/standardized_stats/cluster_binary_F_MNI.nii'], 'file') ~= 2
                  disp(['little output, probably the omnibus contrast from SPM/FAST did not return any voxels for BMMR! subject ' subject]);
                  no_wo_output  = no_wo_output + 1;
                  continue
               end

               cd([path_subjects '/' subject '/standardized_stats']);
               load('power_spectra_one_subject.mat');
               power_spectra = power_spectra + power_spectra_one_subject;

            end

            cd(path_subjects);

            %-average power spectra over subjects
            power_spectra = power_spectra / (no_subjects-no_wo_output);

            system('mv power_spectra.mat power_spectra_old.mat');

            save('power_spectra.mat', 'power_spectra');

         end

      end

   end

end

cd(path_manage)
