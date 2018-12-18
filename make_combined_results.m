

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Combining results from single runs to 'mat' files. Necessary for making figures.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       September 2018 - November 2018
%%%%   Adapted from:  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/make_combined_results_from_single_runs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage         = fgetl(fopen('path_manage.txt'));
path_scratch        = fgetl(fopen('path_scratch.txt'));
path_output         = [path_scratch '/analysis_output_'];
studies_parameters  = readtable([path_manage '/studies_parameters.txt']);
studies             = studies_parameters.study;
packages            = cellstr(['AFNI'; 'FSL '; 'SPM ']);
exper_designs       = cellstr(['boxcar12'; 'boxcar16'; 'boxcar20'; 'event1  '; 'event2  ']);
p                   = parpool(12);
combined_perc       = (-1)*ones(length(studies), length(packages), length(exper_designs), 8, 2, 621);
dims                = size(combined_perc);
pos_rate            = NaN(dims([1 2 3 4 5]));
avg_perc            = NaN(dims([1 2 3 4 5]));
group_perc          = NaN(dims([1 2 3 4 5]));
range_studies       = 1:length(studies);
range_packages      = 1:length(packages);
range_exper_designs = 1:length(exper_designs);


cd(path_manage);
addpath(genpath([path_manage '/matlab_extra_functions']));
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');


%%%%%%%%%-combining 1st level results-%%%%%%%%%%%%%%%%


for study_id   = range_studies

   study       = studies_parameters.study{study_id};
   abbr        = studies_parameters.abbr{study_id};
   task        = studies_parameters.task{study_id};
   no_subjects = studies_parameters.n(study_id);

   for package_id = range_packages

      package     = packages{package_id};

      disp([study ' ' package]);

      if strcmp(package, 'AFNI')
         HRF_models = ['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'tent     '; 'csplin   '];
      elseif strcmp(package, 'FSL')
         HRF_models = ['gamma2   '; 'gamma2_T '; 'gamma    '; 'gamma_T  '; 'FIR      '];
      elseif strcmp(package, 'SPM')
         HRF_models = ['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 '];
      end

      HRF_models    = cellstr(HRF_models);

      parfor subject_id       = 1:no_subjects

         subject              = ['sub-' abbr repmat('0', 1, 4-length(num2str(subject_id))) num2str(subject_id)];
         combined_perc_parfor = (-1)*ones(length(exper_designs), 8, 2);

         for exper_design_id  = range_exper_designs

            exper_design      = exper_designs{exper_design_id};

            for HRF_model_id  = 1:length(HRF_models)

               if (HRF_model_id > 5) && ((exper_design_id ~= 4) || (~strcmp(study, 'CamCAN_sensorimotor')))
                  continue
               end

               HRF_model      = HRF_models{HRF_model_id};

               cd([path_output study '/' package '/exper_design_' exper_design '/HRF_' HRF_model]);

               if exist([subject '/standardized_stats/cluster_binary_F_MNI.nii'], 'file') == 2

                  cd([subject '/standardized_stats']);

                  mask_MNI = niftiread('mask_MNI.nii');

                  if sum(mask_MNI(:)) > 0

                     %-t-test on 1st HRF covariate (canonical function only)
                     if strcmp(HRF_model, 'gamma2') || strcmp(HRF_model, 'gamma2_T') || strcmp(HRF_model, 'gamma2_TD')
                        cluster_binary_t_MNI                                   = niftiread('cluster_binary_t_MNI.nii');
                        combined_perc_parfor(exper_design_id, HRF_model_id, 1) = sum(cluster_binary_t_MNI(:)) / sum(mask_MNI(:));
                     end

                     %-F-test on all HRF covariates
                     cluster_binary_F_MNI                                   = niftiread('cluster_binary_F_MNI.nii');
                     combined_perc_parfor(exper_design_id, HRF_model_id, 2) = sum(cluster_binary_F_MNI(:)) / sum(mask_MNI(:));

                  else

                     disp(['attention! ' pwd]);

                  end

               else

                  disp(['attention! ' pwd]);

               end

            end

         end

         combined_perc(study_id, package_id, :, :, :, subject_id) = combined_perc_parfor;

      end

   end

end

for i1 = 1:dims(1)
   for i2 = 1:dims(2)
      for i3 = 1:dims(3)
         for i4 = 1:dims(4)
            for i5 = 1:dims(5)

               over_sub = combined_perc(i1, i2, i3, i4, i5, :);

               %-(-0.5) chosen for numerical reasons; in fact, we only want to distinguish >=0 from <0
               if sum(over_sub > -0.5) > 0

                  pos_rate(i1, i2, i3, i4, i5) = sum(over_sub > 0) / sum(over_sub > -0.5);
                  %-sum(over_sub < -0.5) considered, as that many times the default (-1) is subtracted; (-1) appears in 'combined_perc' for non-subjects
                  avg_perc(i1, i2, i3, i4, i5) = (sum(over_sub) + sum(over_sub < -0.5)) / sum(over_sub > -0.5);

               end

            end
         end
      end
   end
end

cd(path_manage);
save('combined_results/combined_perc', 'combined_perc');
save('combined_results/pos_rate',      'pos_rate');
save('combined_results/avg_perc',      'avg_perc');

%-saving 'mat' files which show spatial distribution of significant clusters

for study_id      = range_studies

   study          = studies_parameters.study{study_id};
   abbr           = studies_parameters.abbr{study_id};
   no_subjects    = studies_parameters.n(study_id);
   subject_1      = ['sub-' abbr '0001'];

   for package_id = range_packages

      package     = packages{package_id};

      if strcmp(package, 'AFNI')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'tent     '; 'csplin   ']);
      elseif strcmp(package, 'FSL')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma    '; 'gamma_T  '; 'FIR      ']);
      elseif strcmp(package, 'SPM')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 ']);
      end

      parfor exper_design_id = range_exper_designs

         exper_design        = exper_designs{exper_design_id};
         %-loading one mask, only to check the size
         data                = load([path_output study '/' package '/exper_design_' exper_designs{1} '/HRF_gamma2/' subject_1 '/standardized_stats/cluster_binary_F_MNI.mat']);
         dims                = size(data.cluster_binary_F_MNI);
         sp_dist_joint_t     = zeros(length(HRF_models)*dims(1), dims(2), dims(3));
         sp_dist_joint_F     = zeros(length(HRF_models)*dims(1), dims(2), dims(3));

         for HRF_model_id    = 1:5

            HRF_model        = HRF_models{HRF_model_id};

            sp_dist_t        = zeros(dims);
            sp_dist_F        = zeros(dims);

            for subject_id   = 1:no_subjects

               subject                       = ['sub-' abbr repmat('0', 1, 4-length(num2str(subject_id))) num2str(subject_id)];
               cluster_binary_t_MNI_location = [path_output study '/' package '/exper_design_' exper_design '/HRF_' HRF_model '/' subject '/standardized_stats/cluster_binary_t_MNI.mat'];
               cluster_binary_F_MNI_location = [path_output study '/' package '/exper_design_' exper_design '/HRF_' HRF_model '/' subject '/standardized_stats/cluster_binary_F_MNI.mat'];

               if strcmp(HRF_model, 'gamma2') || strcmp(HRF_model, 'gamma2_T') || strcmp(HRF_model, 'gamma2_TD')
                  if exist(cluster_binary_t_MNI_location, 'file') == 2
                     data_t    = load(cluster_binary_t_MNI_location);
                     sp_dist_t = sp_dist_t + data_t.cluster_binary_t_MNI;
                  else
                     disp(['problems with sp_dists_t for ' cluster_binary_t_MNI_location]);
                  end
               end

               if exist(cluster_binary_F_MNI_location, 'file') == 2
                  data_F       = load(cluster_binary_F_MNI_location);
                  sp_dist_F    = sp_dist_F + data_F.cluster_binary_F_MNI;
               else
                  disp(['problems with sp_dists_F for ' cluster_binary_F_MNI_location]);
               end

            end

            if strcmp(HRF_model, 'gamma2') || strcmp(HRF_model, 'gamma2_T') || strcmp(HRF_model, 'gamma2_TD')
               sp_dist_t = rot90(sp_dist_t, 2);
               sp_dist_joint_t((HRF_model_id-1)*dims(1)+1:HRF_model_id*dims(1), :, :) = 100 * sp_dist_t / no_subjects;
            end

            sp_dist_F = rot90(sp_dist_F, 2);
            sp_dist_joint_F((HRF_model_id-1)*dims(1)+1:HRF_model_id*dims(1), :, :) = 100 * sp_dist_F / no_subjects;

         end

         %-trick to save within parfor; parsave2 specifies that the variable is called 'sp_dist_joint'!
         parsave2(['combined_results/sp_dist_joint_t_' study '_' package '_' exper_design '.mat'], sp_dist_joint_t);
         parsave2(['combined_results/sp_dist_joint_F_' study '_' package '_' exper_design '.mat'], sp_dist_joint_F);

      end

   end

end


%%%%%%%%%-combining 2nd level results-%%%%%%%%%%%%%%%%

%-performed separately from single subject analyses, as for single subject analyses, 'parfor' was used before 'exper_designs' and 'HRF_models' loops


for study_id      = range_studies

   study          = studies{study_id};

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

         for HRF_model_id = 1:5

            HRF_model     = HRF_models{HRF_model_id};

            cd([path_output study '/' package '/exper_design_' exper_design '/HRF_' HRF_model]);

            %-t-test on 1st HRF covariate (canonical function only)
            if exist('group_analysis_t_test/indices.txt', 'file') == 2
               indices = textread('group_analysis_t_test/indices.txt');
               mask    = niftiread('mask.nii');
               group_perc(study_id, package_id, exper_design_id, HRF_model_id, 1) = length(indices) / sum(mask(:) == 1);
            end

            %-F-test on all HRF covariates
            if exist('group_analysis_F_test/indices.txt', 'file') == 2
               indices = textread('group_analysis_F_test/indices.txt');
               mask    = niftiread('mask.nii');
               group_perc(study_id, package_id, exper_design_id, HRF_model_id, 2) = length(indices) / sum(mask(:) == 1);
            else
               disp(pwd);
            end


         end

      end

   end

end

cd(path_manage);

save('combined_results/group_perc', 'group_perc');
