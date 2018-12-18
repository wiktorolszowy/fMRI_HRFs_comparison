
%-Wiktor Olszowy

path_manage        = fgetl(fopen('path_manage.txt'));
path_scratch       = fgetl(fopen('path_scratch.txt'));
studies_parameters = readtable([path_manage '/studies_parameters.txt']);
studies            = studies_parameters.study;
packages           = cellstr(['AFNI'; 'FSL '; 'SPM ']);
exper_designs      = cellstr(['boxcar12'; 'boxcar16'; 'boxcar20'; 'event1  '; 'event2  ']);


%-make folders
for study_id = 1:length(studies)

   study     = studies{study_id};
   cd(path_scratch);
   system(['mkdir analysis_output_' study]);

   for package_id = 1:length(packages)

      cd([path_scratch '/analysis_output_' study]);
      package     = packages{package_id};
      system(['mkdir ' package]);
      system('mkdir together_masks');
      path_output = [path_scratch '/analysis_output_' study '/' package];

      if strcmp(package, 'AFNI')
         HRF_models  = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'tent     '; 'csplin   ']);
      elseif strcmp(package, 'FSL')
         HRF_models  = cellstr(['gamma2   '; 'gamma2_T '; 'gamma    '; 'gamma_T  '; 'FIR      ']);
      elseif strcmp(package, 'SPM')
         HRF_models  = cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 ']);
      end

      cd(path_output)
      if strcmp(package, 'FSL')
         system('mkdir preproc_feats_designs');
         system('mkdir preproc_feats');
      elseif strcmp(package, 'SPM')
         system('mkdir preproc');
      end

      for exper_design_id = 1:length(exper_designs)

         exper_design     = exper_designs{exper_design_id};

         cd(path_output)
         system(['mkdir exper_design_', exper_design]);
         path_output_exper_design = [path_output '/exper_design_' exper_design];
         cd(path_output_exper_design);

         for HRF_model_id = 1:length(HRF_models)

            HRF_model     = HRF_models{HRF_model_id};
            system(['mkdir HRF_' HRF_model]);

         end

      end

   end

end

cd(path_manage)
