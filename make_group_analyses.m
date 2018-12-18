

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Group analyses for CamCAN data.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       September 2018 - November 2018
%%%%   Adapted from:  https://github.com/wanderine/ParametricMultisubjectfMRI/blob/master/SPM/run_random_group_analyses_onesamplettest_1.m
%%%%                  http://www.fil.ion.ucl.ac.uk/spm/data/face_rfx/
%%%%                  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/make_group_analyses_random_effects.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage        = fgetl(fopen('path_manage.txt'));
path_scratch       = fgetl(fopen('path_scratch.txt'));
studies_parameters = readtable([path_manage '/studies_parameters.txt']);
studies            = studies_parameters.study;
packages           = cellstr(['AFNI'; 'FSL '; 'SPM ']);
exper_designs      = cellstr(['boxcar12'; 'boxcar16'; 'boxcar20'; 'event1  '; 'event2  ']);

addpath(genpath([path_manage '/matlab_extra_functions']));
addpath('/applications/spm/spm12_7219');
addpath('/applications/AFNI/AFNI_18.0.11');

%-following Eklund ea 2016, I only consider sample size 20
no_subjects        = 20;

spm('Defaults', 'fMRI');
spm_jobman('initcfg');


for study_id = 1:length(studies)

   study     = studies{study_id};

   for package_id = 1:length(packages)

      package     = packages{package_id};

      if strcmp(package, 'AFNI')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'tent     '; 'csplin   ']);
      elseif strcmp(package, 'FSL')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma    '; 'gamma_T  '; 'FIR      ']);
      elseif strcmp(package, 'SPM')
         HRF_models = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 ']);
      end

      for exper_design_id = 1:length(exper_designs)

         exper_design     = exper_designs{exper_design_id};

         for HRF_model_id = 1:length(HRF_models)

            HRF_model     = HRF_models{HRF_model_id};
            path_output   = [path_scratch '/analysis_output_' study '/' package '/exper_design_' exper_design '/HRF_' HRF_model];

            %-weird numerical problems
            if strcmp(study, 'BMMR_checkerboard') && strcmp(package, 'SPM')  && strcmp(exper_design, 'boxcar12') && strcmp(HRF_model, 'Fourier')
               continue
            end

            %-no voxels with initial zstat1>3.1 -> numerical problems
            if strcmp(study, 'CRIC_checkerboard') && strcmp(package, 'AFNI') && strcmp(exper_design, 'boxcar20') && strcmp(HRF_model, 'csplin')
               continue
            end
            if strcmp(study, 'CRIC_checkerboard') && strcmp(package, 'AFNI') && strcmp(exper_design, 'event1')   && strcmp(HRF_model, 'tent')
               continue
            end
            if strcmp(study, 'CRIC_checkerboard') && strcmp(package, 'AFNI') && strcmp(exper_design, 'event1')   && strcmp(HRF_model, 'csplin')
               continue
            end
            if strcmp(study, 'CRIC_checkerboard') && strcmp(package, 'SPM')  && strcmp(exper_design, 'event1')   && strcmp(HRF_model, 'gamma2_T')
               continue
            end
            if strcmp(study, 'BMMR_checkerboard') && strcmp(package, 'SPM')  && strcmp(exper_design, 'boxcar16') && strcmp(HRF_model, 'gamma_x2')
               continue
            end

            %-different numbers of regressors/parameters of interest for the different HRF models
            if strcmp(HRF_model, 'gamma') || strcmp(HRF_model, 'gamma2')
               no_par = 1;
            elseif strcmp(HRF_model, 'gamma_T') || strcmp(HRF_model, 'gamma2_T')
               no_par = 2;
            elseif strcmp(HRF_model, 'gamma2_TD')
               no_par = 3;
            elseif strcmp(HRF_model, 'Fourier')
               no_par = 11;
            elseif strcmp(HRF_model, 'FIR') || strcmp(HRF_model, 'tent') || strcmp(HRF_model, 'csplin')
               no_par = 6;
            elseif strcmp(HRF_model, 'gamma_x1')
               no_par = 1;
            elseif strcmp(HRF_model, 'gamma_x2')
               no_par = 2;
            elseif strcmp(HRF_model, 'gamma_x3')
               no_par = 3;
            end

            %-when the 'group_analysis' folder is not empty, problems occur
            %-when repeating a run, slight differences in results, probably some randomization involved (?)

            %%%%%%%%%%% F-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            cd(path_output);

            system('rm -rf group_analysis_F_test');
            system('mkdir group_analysis_F_test');
            clear jobs;
            SPM_mat_location = cellstr(fullfile(path_output, 'group_analysis_F_test', 'SPM.mat'));
            
            jobs{1}.stats{1}.factorial_design.dir                            = cellstr(fullfile(path_output, 'group_analysis_F_test'));
            %-it seems there are errors when the name is not specified
            jobs{1}.stats{1}.factorial_design.des.fd.fact.name               = 'Basis';
            jobs{1}.stats{1}.factorial_design.des.fd.fact.levels             = no_par;
            jobs{1}.stats{1}.factorial_design.des.fd.fact.dept               = 1;
            for par_id = 1:no_par
               jobs{1}.stats{1}.factorial_design.des.fd.icell(par_id).levels = par_id;
               coef_maps                                                     = cellstr(spm_select('FPListRec', fullfile(path_output), ['^coef_MNI_masked_' num2str(par_id) '.nii']));
               jobs{1}.stats{1}.factorial_design.des.fd.icell(par_id).scans  = coef_maps(1:no_subjects);
            end

            %-creating brain mask over all subjects (up to 'no_subjects')
            AFNI_mask_command    = ['3dMean -prefix ' path_output '/mask.nii -mask_inter'];
            for subject_id       = 1:no_subjects
               AFNI_mask_command = [AFNI_mask_command ' ' coef_maps{subject_id}];
            end
            system('rm mask.nii');
            system(AFNI_mask_command);

            jobs{1}.stats{1}.factorial_design.masking.em                     = {[path_output '/mask.nii']};
            jobs{1}.stats{2}.fmri_est.spmmat                                 = SPM_mat_location;
            jobs{1}.stats{3}.con.spmmat                                      = SPM_mat_location;
            jobs{1}.stats{3}.con.consess{1}.fcon.name                        = 'F-contrast';
            jobs{1}.stats{3}.con.consess{1}.fcon.weights                     = eye(no_par);
            jobs{1}.stats{4}.results.spmmat                                  = SPM_mat_location;
            jobs{1}.stats{4}.results.conspec.contrasts                       = Inf;
            jobs{1}.stats{4}.results.conspec.threshdesc                      = 'FWE';
            
            spm_jobman('run', jobs);
            
            %-get F-map
            V           = spm_vol([path_output '/group_analysis_F_test/spmF_0001.nii']);
            [tmap,aa]   = spm_read_vols(V);
            
            %-calculate cluster extent threshold
            [k,Pc]      = corrclusth(SPM, 0.001, 0.05, 1:100000);
            df          = [no_par SPM.xX.erdf];
            u           = spm_u(0.001, df, 'F');
            indices     = find(tmap>u);
            
            %-get the size of the largest cluster
            max_cluster = max_extent(tmap, indices);
            
            if max_cluster >= k
               indices     = find(tmap>u);
            else
               indices     = '';
            end
            fid            = fopen('indices.txt', 'wt');
            fprintf(fid, '%8.0f\n', indices);

            %%%%%%%%%%% t-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if strcmp(HRF_model, 'gamma2') || strcmp(HRF_model, 'gamma2_T') || strcmp(HRF_model, 'gamma2_TD')

               cd(path_output);

               system('rm -rf group_analysis_t_test');
               system('mkdir group_analysis_t_test');
               clear jobs;
               SPM_mat_location = cellstr(fullfile(path_output, 'group_analysis_t_test', 'SPM.mat'));

               coef_maps        = cellstr(spm_select('FPListRec', fullfile(path_output), '^coef_MNI_masked_1.nii'));
               
               jobs{1}.stats{1}.factorial_design.dir            = cellstr(fullfile(path_output, 'group_analysis_t_test'));
               jobs{1}.stats{1}.factorial_design.des.t1.scans   = coef_maps(1:no_subjects);
               jobs{1}.stats{1}.factorial_design.cov            = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
               %-using the explicit mask (generated for F-test)
               jobs{1}.stats{1}.factorial_design.masking.em     = {[path_output '/mask.nii']};
               jobs{1}.stats{2}.fmri_est.spmmat                 = SPM_mat_location;

               jobs{1}.stats{3}.con.spmmat                      = SPM_mat_location;
               jobs{1}.stats{3}.con.consess{1}.tcon             = struct('name', 't-contrast', 'convec', 1, 'sessrep', 'none');

               jobs{1}.stats{4}.results.spmmat                  = SPM_mat_location;
               jobs{1}.stats{4}.results.conspec.contrasts       = Inf;
               jobs{1}.stats{4}.results.conspec.threshdesc      = 'FWE';
               
               spm_jobman('run', jobs);

               %-get t-map
               V           = spm_vol([path_output '/group_analysis_t_test/spmT_0001.nii']);
               [tmap,aa]   = spm_read_vols(V);

               %-calculate cluster extent threshold
               [k,Pc]      = corrclusth(SPM, 0.001, 0.05, 1:100000);
               df          = [1 SPM.xX.erdf];
               u           = spm_u(0.001, df, 'T');
               indices     = find(tmap>u);
               
               %-get the size of the largest cluster
               max_cluster = max_extent(tmap, indices);
               
               if max_cluster >= k
                  indices     = find(tmap>u);
               else
                  indices     = '';
               end
               fid            = fopen('indices.txt', 'wt');
               fprintf(fid, '%8.0f\n', indices);

            end

            %-saving 3D maps showing where the significant clusters are located

            cd([path_output '/group_analysis_F_test']);
            mask                    = niftiread('mask.nii');
            cluster_binary          = zeros(size(mask));
            indices                 = fscanf(fopen('indices.txt', 'r'), '%d');
            cluster_binary(indices) = 1;
            niftiwrite(cluster_binary, 'cluster_binary.nii');

            if exist([path_output '/group_analysis_t_test'], 'dir') == 7
               cd([path_output '/group_analysis_t_test']);
               mask                    = niftiread('mask.nii');
               cluster_binary          = zeros(size(mask));
               indices                 = fscanf(fopen('indices.txt', 'r'), '%d');
               cluster_binary(indices) = 1;
               niftiwrite(cluster_binary, 'cluster_binary.nii');
            end

         end

      end

   end

end

cd(path_manage)
