

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Figures comparing results for different HRF models.
%%%%   Written by:  Wiktor Olszowy, University of Cambridge
%%%%   Contact:     wo222@cam.ac.uk
%%%%   Created:     September 2018 - November 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


paper                = 'HRF_comp';
path_manage          = fgetl(fopen('path_manage.txt'));
path_scratch         = fgetl(fopen('path_scratch.txt'));
path_output          = [path_scratch '/analysis_output_'];
studies_parameters   = readtable([path_manage '/studies_parameters.txt']);
studies              = studies_parameters.study;
packages             = cellstr(['AFNI'; 'FSL '; 'SPM ']);
exper_designs        = cellstr(['boxcar12'; 'boxcar16'; 'boxcar20'; 'event1  '; 'event2  ']);
clims                = [80 80 80 30];
exper_designs_exp_id = [4 2 3 3 1];
freq_studies_exp_id  = [100 1/32 1/40 1/40 1/24];
colors               = [0 1 0; 0.96 0.47 0.13; 1 0 1; 0 1 1; 0 0 1];
range_packages       = 1:length(packages);
range_studies        = 1:length(studies);
range_exper_designs  = 1:length(exper_designs);
slices               = [25 37 54 71];
clims_studies        = [70 20 20 20 70];
dims_MNI             = [91 109 91];

cd(path_manage);
addpath(genpath([path_manage '/matlab_extra_functions']));
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');

load('combined_results/pos_rate.mat');
load('combined_results/avg_perc.mat');
load('combined_results/group_perc.mat');

studies_labels = studies;
for study_id = range_studies
   study       = studies{study_id};
   TR          = studies_parameters.TR(study_id);
   study_label = strrep(study, '_', ' ');
   study_label = strrep(study_label, '1400', '');
   study_label = strrep(study_label, '645',  '');
   study_label = strrep(study_label, ' release 3',  '');
   studies_labels{study_id} = [study_label ' (TR=' num2str(TR) 's)'];
end

HRF_models        = [{cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'tent     '; 'csplin   '])},                                           %-AFNI
                     {cellstr(['gamma2   '; 'gamma2_T '; 'gamma    '; 'gamma_T  '; 'FIR      '])},                                           %-FSL
                     {cellstr(['gamma2   '; 'gamma2_T '; 'gamma2_TD'; 'Fourier  '; 'FIR      '; 'gamma_x1 '; 'gamma_x2 '; 'gamma_x3 '])}];   %-SPM
HRF_models_labels = [{cellstr(['gam2     '; 'gam2+T   '; 'gam2+TD  '; 'tent     '; 'csplin   '])},                                           %-AFNI
                     {cellstr(['gam2     '; 'gam2+T   '; 'gam1     '; 'gam1+T   '; 'FIR      '])},                                           %-FSL
                     {cellstr(['gam2     '; 'gam2+T   '; 'gam2+TD  '; 'Fourier  '; 'FIR      '; 'gamma x 1'; 'gamma x 2'; 'gamma x 3'])}];   %-SPM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st LEVEL: POWER SPECTRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrows    = 6;
fig_size = [0 0 450 1.1*695];
study_id_plot     = 0;
figure('rend', 'painters', 'pos', fig_size, 'Visible', 'off');
for study_id      = range_studies
   study          = studies{study_id};
   study_label    = studies_labels{study_id};
   %-omitting CRIC checkerboard as group brain mask without V1
   no_subjects    = studies_parameters.n(study_id);
   TR             = studies_parameters.TR(study_id);
   exper_design   = exper_designs{exper_designs_exp_id(study_id)};
   power_spectra_all = (-1)*ones(5, 5, 257);
   for package_id = range_packages
      for HRF_model_id = 1:5
         power_spectra = load([path_output study '/' packages{package_id} '/exper_design_' exper_design '/HRF_' HRF_models{package_id}{HRF_model_id} '/power_spectra.mat']);
         power_spectra_all(package_id, HRF_model_id, 1:257) = power_spectra.power_spectra(1:257);
      end
   end
   for package_id = range_packages
      package       = packages{package_id};
      study_id_plot = study_id_plot + 1;
      subplot(nrows, 3, study_id_plot);
      f = linspace(0, 0.5/TR, 257);
      axis off;
      %-[left bottom width height]
      gca_pos = get(gca, 'Position');
      %-lowering vertical spacing between subplots
      if study_id > 1
         gca_pos(2) = gca_pos(2) + (study_id-1)*0.066;
      end
      %-lowering horizontal spacing between subplots
      if package_id > 1
         gca_pos(1) = gca_pos(1) - (package_id-1)*0.08;
      end
      ax1     = axes('Position', gca_pos, 'Visible', 'off');
      gca_pos = [gca_pos(1)+0.05  gca_pos(2)+0.05  0.92*gca_pos(3)  gca_pos(4)-0.05];
      ax2     = axes('Position', gca_pos, 'Visible', 'off');
      h1  = plot(f, reshape(power_spectra_all(package_id, 1, :), [1 257]), '-.', 'color', colors(1, :)); hold on;
      h2  = plot(f, reshape(power_spectra_all(package_id, 2, :), [1 257]), '-.', 'color', colors(2, :)); hold on;
      h3  = plot(f, reshape(power_spectra_all(package_id, 3, :), [1 257]), '-.', 'color', colors(3, :)); hold on;
      h4  = plot(f, reshape(power_spectra_all(package_id, 4, :), [1 257]), '-.', 'color', colors(4, :)); hold on;
      h5  = plot(f, reshape(power_spectra_all(package_id, 5, :), [1 257]), '-.', 'color', colors(5, :)); hold on;
            plot(freq_studies_exp_id(study_id), 0,                         'kx', 'MarkerSize', 4); hold on;
      h10 = plot(freq_studies_exp_id(study_id), max(power_spectra_all(:)), 'kx', 'MarkerSize', 4); hold on;
      h80 = plot([0 0.5/TR], [1 1], 'k');                                      hold on;
      set([h1 h2 h3 h4 h5 h80], 'LineWidth', 0.8);
      if study_id == 5
         hx = xlabel('Frequency [Hz]');
      else
         hx = '';
      end
      if package_id == 1
         hy = ylabel('Power spectra', 'Units', 'normalized');
      else
         hy = '';
      end
      htitle = title('');
      set([hx hy htitle], 'FontSize', 5);
      if study_id == 1 && package_id == 2
         htitle = title({package; ' '; study_label});
      elseif study_id == 1
         htitle = title({package; ' '; ' '});
      elseif package_id == 2
         htitle = title({' '; study_label});
      end
      xlim([0 0.5/TR]);
      ylim([0 max(power_spectra_all(:))]);
      %-controlling distance between y axis title and y axis
      hy_pos    = get(hy, 'Position');
      hy_pos(1) = -0.11;
      set(hy, 'Position', hy_pos);
      axes(ax2);
      x_ax_marks = linspace(0, 0.5/TR, 6);
      if package_id == 1
         set(gca, 'XTick', x_ax_marks);
         set(gca, 'XTickLabel', round(x_ax_marks, 2), 'FontSize', 3.5);
      else
         %-omitting zeros, as there is very little horizontal spacing between the subplots
         set(gca, 'XTick', x_ax_marks(2:length(x_ax_marks)));
         set(gca, 'XTickLabel', round(x_ax_marks(2:length(x_ax_marks)), 2), 'FontSize', 3.5);         
         set(gca, 'YTickLabel', []);
      end
      %-trick to specify the fontsize again, so that fontsize 3.5 for 'XTickLabel' doesn't affect the y-axis label
      set([hx hy htitle], 'FontSize', 5);
   end
   fig_ref = subplot(nrows, 3, nrows*3-2);
   plot(1);
   set(gca, 'visible', 'off');
   %-buffer values: first value down -> text to the left; second value down -> text down
   legendflex([h1 h2 h3 h4 h5], {HRF_models_labels{1}{1}, HRF_models_labels{1}{2}, HRF_models_labels{1}{3}, HRF_models_labels{1}{4}, HRF_models_labels{1}{5}}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-59  275], 'box', 'off', 'FontSize', 4, 'nrow', 5);
   fig_ref = subplot(nrows, 3, nrows*3-1);
   plot(1);
   set(gca, 'visible', 'off');
   legendflex([h1 h2 h3 h4 h5], {HRF_models_labels{2}{1}, HRF_models_labels{2}{2}, HRF_models_labels{2}{3}, HRF_models_labels{2}{4}, HRF_models_labels{2}{5}}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-94 275], 'box', 'off', 'FontSize', 4, 'nrow', 5);
   fig_ref = subplot(nrows, 3, nrows*3);
   plot(1);
   set(gca, 'visible', 'off');
   legendflex([h1 h2 h3 h4 h5], {HRF_models_labels{3}{1}, HRF_models_labels{3}{2}, HRF_models_labels{3}{3}, HRF_models_labels{3}{4}, HRF_models_labels{3}{5}}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-130 275], 'box', 'off', 'FontSize', 4, 'nrow', 5);
end
figname = [paper '_power'];
print_to_svg_to_pdf(figname, path_manage);


%%%%%%%%%%%%%%%%%%%%%%%% 1st LEVEL: QQ-PLOTS: 1ST SUBJECT IN EACH DATASET ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
nrows    = 6;
fig_size = [0 0 450 1.1*695];
study_id_plot     = 0;
figure('rend', 'painters', 'pos', fig_size, 'Visible', 'off');
for study_id      = range_studies
   study          = studies{study_id};
   study_label    = studies_labels{study_id};
   %-omitting CRIC checkerboard as group brain mask without V1
   no_subjects    = studies_parameters.n(study_id);
   TR             = studies_parameters.TR(study_id);
   abbr           = studies_parameters.abbr{study_id};
   subject_1      = ['sub-' abbr '0001'];
   exper_design   = exper_designs{exper_designs_exp_id(study_id)};
   %-0.3: upper boundary for the proportion of brain voxels; 91x109x91: MNI dimensions; 261: max number of time points
   res4d_all      = (-1)*ones(3, 5, round(0.3*91*109*91*261));
   for package_id = range_packages
      for HRF_model_id  = 1:5
         res4d_location = [path_output study '/' packages{package_id} '/exper_design_' exper_design '/HRF_' HRF_models{package_id}{HRF_model_id} '/', subject_1 '/standardized_stats'];
         system(['gunzip ' res4d_location '/res4d_masked.nii.gz']);
         res4d   = niftiread([res4d_location '/res4d_masked.nii']);
         system(['gzip ' res4d_location '/res4d_masked.nii']);
         %-normalize the signal in each voxel separately
         res4d   = zscore(res4d, [], 4);
         res4d   = reshape(res4d, 1, []);
         res4d   = res4d(res4d~=0);
         res4d   = res4d(~isnan(res4d));
         res4d   = zscore(res4d);
         res4d_all(package_id, HRF_model_id, 1:length(res4d)) = res4d;
      end
   end
   %-control random number generation
   rng(1);
   %-normal theoretical quantiles
   norm_th_q = sort(normrnd(0, 1, [1 length(res4d)]));
   for package_id   = range_packages
      package       = packages{package_id};
      study_id_plot = study_id_plot + 1;
      subplot(nrows, 3, study_id_plot);
      axis off;
      %-[left bottom width height]
      gca_pos = get(gca, 'Position');
      %-lowering vertical spacing between subplots
      if study_id > 1
         gca_pos(2) = gca_pos(2) + (study_id-1)*0.066;
      end
      %-lowering horizontal spacing between subplots
      if package_id > 1
         gca_pos(1) = gca_pos(1) - (package_id-1)*0.08;
      end
      ax1     = axes('Position', gca_pos, 'Visible', 'off');
      gca_pos = [gca_pos(1)+0.05  gca_pos(2)+0.05  0.92*gca_pos(3)  gca_pos(4)-0.05];
      ax2     = axes('Position', gca_pos, 'Visible', 'off');
      h1      = plot(norm_th_q, sort(reshape(res4d_all(package_id, 1, 1:length(res4d)), [1 length(res4d)])), '-.', 'color', colors(1, :)); hold on;
      h2      = plot(norm_th_q, sort(reshape(res4d_all(package_id, 2, 1:length(res4d)), [1 length(res4d)])), '-.', 'color', colors(2, :)); hold on;
      h3      = plot(norm_th_q, sort(reshape(res4d_all(package_id, 3, 1:length(res4d)), [1 length(res4d)])), '-.', 'color', colors(3, :)); hold on;
      h4      = plot(norm_th_q, sort(reshape(res4d_all(package_id, 4, 1:length(res4d)), [1 length(res4d)])), '-.', 'color', colors(4, :)); hold on;
      h5      = plot(norm_th_q, sort(reshape(res4d_all(package_id, 5, 1:length(res4d)), [1 length(res4d)])), '-.', 'color', colors(5, :)); hold on;
      %-ideal fit with normal distribution
      h10     = plot([-1000 1000], [-1000 1000],  'k'); hold on;
      set([h1 h2 h3 h4 h5 h10], 'LineWidth', 0.8);
      if study_id == 5
         hx = xlabel('Normal theoretical quantiles');
      else
         hx = '';
      end
      if package_id == 1
         hy = ylabel('Data quantiles', 'Units', 'normalized');
      else
         hy = '';
      end
      set([hx hy htitle], 'FontSize', 5);
      if study_id == 1 && package_id == 2
         htitle = title({package; ' '; study_label});
      elseif study_id == 1
         htitle = title({package; ' '; ' '});
      elseif package_id == 2
         htitle = title({' '; study_label});
      end
      xlim([min(norm_th_q) max(norm_th_q)]);
      ylim([min(res4d_all(:)) max(res4d_all(:))]);
      %-controlling distance between y axis title and y axis
      hy_pos    = get(hy, 'Position');
      hy_pos(1) = -0.11;
      set(hy, 'Position', hy_pos);
      axes(ax2);
      a = get(gca, 'XTickLabel');
      set(gca, 'XTickLabel', a, 'FontSize', 3.5);
      if package_id > 1
         set(gca, 'YTickLabel', []);
      end
      %-trick to specify the fontsize again, so that fontsize 3.5 for 'XTickLabel' doesn't affect the y-axis label
      set([hx hy htitle], 'FontSize', 5);
   end
   fig_ref = subplot(nrows, 3, nrows*3-2);
   plot(1);
   set(gca, 'visible', 'off');
   %-buffer values: first value down -> text to the left; second value down -> text down
   legendflex([h1 h2 h3 h4 h5], {HRF_models_labels{1}{1}, HRF_models_labels{1}{2}, HRF_models_labels{1}{3}, HRF_models_labels{1}{4}, HRF_models_labels{1}{5}}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-59  275], 'box', 'off', 'FontSize', 4, 'nrow', 5);
   fig_ref = subplot(nrows, 3, nrows*3-1);
   plot(1);
   set(gca, 'visible', 'off');
   legendflex([h1 h2 h3 h4 h5], {HRF_models_labels{2}{1}, HRF_models_labels{2}{2}, HRF_models_labels{2}{3}, HRF_models_labels{2}{4}, HRF_models_labels{2}{5}}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-94 275], 'box', 'off', 'FontSize', 4, 'nrow', 5);
   fig_ref = subplot(nrows, 3, nrows*3);
   plot(1);
   set(gca, 'visible', 'off');
   legendflex([h1 h2 h3 h4 h5], {HRF_models_labels{3}{1}, HRF_models_labels{3}{2}, HRF_models_labels{3}{3}, HRF_models_labels{3}{4}, HRF_models_labels{3}{5}}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-130 275], 'box', 'off', 'FontSize', 4, 'nrow', 5);
end
figname = [paper '_qqplots_1st_subject'];
print_to_svg_to_pdf(figname, path_manage);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st+2nd LEVELS: PERCENTAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-level 1: single subject
%-level 2: group
parfor level_type = [1 2]
   %-test 1: t-test
   %-test 2: F-test
   for test_type = [1 2]
      if test_type == 1
         test = 't';
      else
         test = 'F';
      end
      if level_type == 1
         level    = '1st';
         results  = avg_perc;
         ylab     = 'Avg. % of sig voxels';
         nrows    = 6;
         fig_size = [0 0 450 1.1*695];
      else
         level    = '2nd';
         results  = group_perc;
         ylab     = '% of sig voxels';
         nrows    = 5;
         fig_size = [0 0 450 695];
      end
      study_id_plot     = 0;
      figure('rend', 'painters', 'pos', fig_size, 'Visible', 'off');
      for study_id      = range_studies
         study          = studies{study_id};
         study_label    = studies_labels{study_id};
         %-omitting CRIC checkerboard as group brain mask without V1
         if (strcmp(level, '2nd')) && (strcmp(study, 'CRIC_checkerboard'))
            continue
         end
         no_subjects    = studies_parameters.n(study_id);
         for package_id = range_packages
            package       = packages{package_id};
            study_id_plot = study_id_plot + 1;
            subplot(nrows, 3, study_id_plot);
            axis off;
            %-[left bottom width height]
            gca_pos = get(gca, 'Position');
            %-lowering vertical spacing between subplots
            if level_type == 1
               if study_id > 1
                  gca_pos(2) = gca_pos(2) + (study_id-1)*0.056;
               end
            elseif study_id > 2
               gca_pos(2) = gca_pos(2) + (study_id-2)*0.056;
            end
            %-lowering horizontal spacing between subplots
            if package_id > 1
               gca_pos(1) = gca_pos(1) - (package_id-1)*0.08;
            end
            ax1     = axes('Position', gca_pos, 'Visible', 'off');
            gca_pos = [gca_pos(1)+0.05  gca_pos(2)+0.05  0.92*gca_pos(3)  gca_pos(4)-0.05];
            ax2     = axes('Position', gca_pos, 'Visible', 'off');
            if test_type == 2
               HRF_models_length = 5;
            elseif strcmp(package, 'AFNI') || strcmp(package, 'SPM')
               HRF_models_length = 3;
            else
               HRF_models_length = 2;
            end
            %-ylim up to above max_y (*1.03), so that there is some space above the highest bar
            max_y = 1.03*max(100*reshape(results(study_id, 1:3, 1:5, 1:5, 1:2), 3*5*5*2, 1));
            for HRF_model_id = 1:HRF_models_length
               x_axis = HRF_model_id + (HRF_model_id-1)*length(exper_designs);
               if test_type == 1
                  y_axis = 100*results(study_id, package_id, :, HRF_model_id, 1);
               else
                  y_axis = 100*results(study_id, package_id, :, HRF_model_id, 2);
               end
               h1     = bar(ax2, x_axis+1, y_axis(1), 'FaceColor', colors(1,:), 'LineWidth', 0.5); hold on;
               h2     = bar(ax2, x_axis+2, y_axis(2), 'FaceColor', colors(2,:), 'LineWidth', 0.5); hold on;
               h3     = bar(ax2, x_axis+3, y_axis(3), 'FaceColor', colors(3,:), 'LineWidth', 0.5); hold on;
               h4     = bar(ax2, x_axis+4, y_axis(4), 'FaceColor', colors(4,:), 'LineWidth', 0.5); hold on;
               h5     = bar(ax2, x_axis+5, y_axis(5), 'FaceColor', colors(5,:), 'LineWidth', 0.5); hold on;
               %-marking the true design
               dist_5 = 5*(max_y/100);
               if test_type == 1
                  fill = (4.5*dist_5/5:dist_5:(y_axis(exper_designs_exp_id(study_id))-4.0*dist_5/5));
               else
                  fill = (3*  dist_5/5:dist_5:(y_axis(exper_designs_exp_id(study_id))-1.5*dist_5/5));
               end
               h16    = plot(transpose(repmat(x_axis+exper_designs_exp_id(study_id), length(fill), 1)), fill, '+'); hold on;
               if test_type == 1
                  set(h16, 'markersize', 2.4, 'Color', 'black');
               else
                  set(h16, 'markersize', 1.5, 'Color', 'black');
               end
            end
            if package_id == 1
               hy = ylabel(ylab, 'Units', 'normalized');
            else
               hy = '';
            end
            if level_type == 1
               if study_id == 1 && package_id == 2
                  htitle = title({package; ' '; study_label});
               elseif study_id == 1
                  htitle = title({package; ' '; ' '});
               elseif package_id == 2
                  htitle = title({' '; study_label});
               end
            else
               if study_id == 1 && package_id == 2
                  htitle = title({package; ' '; study_label; ' '});
               elseif study_id == 1
                  htitle = title({package; ' '; ' '; ' '});
               elseif package_id == 2
                  htitle = title({' '; study_label; ' '});
               end
            end
            if test_type == 1
               xlim([0.8 (5+1)*3+1+0.2]);
            else
               xlim([0.8 (5+1)*5+1+0.2]);
            end
            %-upper limit for ylim cannot be equal lower limit
            if max_y > 0
               ylim([0 max_y]);
            else
               ylim([0 2]);
            end
            %-controlling distance between y axis title and y axis
            hy_pos    = get(hy, 'Position');
            hy_pos(1) = -0.11;
            set(hy, 'Position', hy_pos);
            axes(ax2);
            set(gca, 'XTick', 1000); %-1000: just some value beyond x-axis
            set(gca, 'XTickLabel', repmat({' '}, 1, 1), 'FontSize', 3);
            if package_id > 1
               set(gca, 'YTickLabel', []);
            end
            set([hy htitle], 'FontSize', 5);
            %-you have to find the two values manually (middle locations of first and last HRF labels)
            if test_type == 1
               %text(ax1, 0.395, 0.2, 'z', 'HorizontalAlignment', 'center', 'fontweight', 'bold', 'FontSize', 4);
               %text(ax1, 0.990, 0.2, 'z', 'HorizontalAlignment', 'center', 'fontweight', 'bold', 'FontSize', 4);
               from = 0.395;
               to   = 0.990;
            else
               %text(ax1, 0.335, 0.2, 'z', 'HorizontalAlignment', 'center', 'fontweight', 'bold', 'FontSize', 4);
               %text(ax1, 1.055, 0.2, 'z', 'HorizontalAlignment', 'center', 'fontweight', 'bold', 'FontSize', 4);
               from = 0.335;
               to   = 1.055;
            end
            for HRF_model_id = 1:HRF_models_length
               %-specifying font size
               if test_type == 1
                  n  = 3;
                  fs = 4.5;
               else
                  n  = 5;
                  fs = 3.25;
               end
               if level_type == 1
                  text(ax1, from+(to-from)*(HRF_model_id-1)/(n-1), 0.36, HRF_models_labels{package_id}{HRF_model_id}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', fs, 'Units', 'normalized');
               else
                  text(ax1, from+(to-from)*(HRF_model_id-1)/(n-1), 0.30, HRF_models_labels{package_id}{HRF_model_id}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', fs, 'Units', 'normalized');
               end
            end   
         end
         fig_ref = subplot(nrows, 3, nrows*3);
         plot(1);
         set(gca, 'visible', 'off');
         %-buffer values: first value down -> text to the left; second value down -> text down
         if level_type == 1
            hlegend = legendflex([h1 h2 h3 h4 h5 h16], {sprintf('boxcar 12s off + 12s on') sprintf('boxcar 16s off + 16s on') sprintf('boxcar 20s off + 20s on') sprintf('CamCAN event-related design') sprintf('~U(3,6)s off + 0.1s on') 'True experimental design'}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-342 250], 'box', 'off', 'FontSize', 4, 'nrow', 2);
         else
            hlegend = legendflex([h1 h2 h3 h4 h5 h16], {sprintf('boxcar 12s off + 12s on') sprintf('boxcar 16s off + 16s on') sprintf('boxcar 20s off + 20s on') sprintf('CamCAN event-related design') sprintf('~U(3,6)s off + 0.1s on') 'True experimental design'}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-342 200], 'box', 'off', 'FontSize', 4, 'nrow', 2);
         end
      end
      figname = [paper '_' level '_level_percentages_' test];
      print_to_svg_to_pdf(figname, path_manage);
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st LEVEL: POSITIVE RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-test 1: t-test
%-test 2: F-test
for test_type   = [1 2]
   if test_type == 1
      test  = 't';
   else
      test  = 'F';
   end
   nrows    = 6;
   fig_size = [0 0 450 1.1*695];
   study_id_plot     = 0;
   figure('rend', 'painters', 'pos', fig_size, 'Visible', 'off');
   for study_id      = range_studies
      study          = studies{study_id};
      study_label    = studies_labels{study_id};
      %-omitting CRIC checkerboard as group brain mask without V1
      no_subjects    = studies_parameters.n(study_id);
      %-needed for the CIs
      sd             = sqrt(0.05*0.95/no_subjects)*1.96*100;
      for package_id = range_packages
         package       = packages{package_id};
         study_id_plot = study_id_plot + 1;
         subplot(nrows, 3, study_id_plot);
         axis off;
         %-[left bottom width height]
         gca_pos = get(gca, 'Position');
         %-lowering vertical spacing between subplots
         if study_id > 1
            gca_pos(2) = gca_pos(2) + (study_id-1)*0.056;
         end
         %-lowering horizontal spacing between subplots
         if package_id > 1
            gca_pos(1) = gca_pos(1) - (package_id-1)*0.08;
         end
         ax1     = axes('Position', gca_pos, 'Visible', 'off');
         gca_pos = [gca_pos(1)+0.05  gca_pos(2)+0.05  0.92*gca_pos(3)  gca_pos(4)-0.05];
         ax2     = axes('Position', gca_pos, 'Visible', 'off');
         if test_type == 2
            HRF_models_length = 5;
         elseif strcmp(package, 'AFNI') || strcmp(package, 'SPM')
            HRF_models_length = 3;
         else
            HRF_models_length = 2;
         end
         %-ylim up to above max_y (*1.03), so that there is some space above the highest bar
         max_y = 1.03*max(100*reshape(pos_rate(study_id, 1:3, 1:5, 1:5, 1:2), 3*5*5*2, 1));
         for HRF_model_id = 1:HRF_models_length
            x_axis = HRF_model_id + (HRF_model_id-1)*length(exper_designs);
            if test_type == 1
               y_axis = 100*pos_rate(study_id, package_id, :, HRF_model_id, 1);
            else
               y_axis = 100*pos_rate(study_id, package_id, :, HRF_model_id, 2);
            end
            h1     = bar(ax2, x_axis+1, y_axis(1), 'FaceColor', colors(1,:), 'LineWidth', 0.5); hold on;
            h2     = bar(ax2, x_axis+2, y_axis(2), 'FaceColor', colors(2,:), 'LineWidth', 0.5); hold on;
            h3     = bar(ax2, x_axis+3, y_axis(3), 'FaceColor', colors(3,:), 'LineWidth', 0.5); hold on;
            h4     = bar(ax2, x_axis+4, y_axis(4), 'FaceColor', colors(4,:), 'LineWidth', 0.5); hold on;
            h5     = bar(ax2, x_axis+5, y_axis(5), 'FaceColor', colors(5,:), 'LineWidth', 0.5); hold on;
            %-marking the true design
            dist_5 = 5*(max_y/100);
            if test_type == 1
               fill = (4.5*dist_5/5:dist_5:(y_axis(exper_designs_exp_id(study_id))-4.0*dist_5/5));
            else
               fill = (3*  dist_5/5:dist_5:(y_axis(exper_designs_exp_id(study_id))-1.5*dist_5/5));
            end
            h16    = plot(transpose(repmat(x_axis+exper_designs_exp_id(study_id), length(fill), 1)), fill, '+'); hold on;
            if test_type == 1
               set(h16, 'markersize', 2.4, 'Color', 'black');
            else
               set(h16, 'markersize', 1.5, 'Color', 'black');
            end
         end
         %-confidence interval (around 5%)
         h26 = plot([0 1000], [5    5],    'k', 'LineWidth', 1.5, 'Color', [0.2 0   0]); hold on;   %-brown
         h27 = plot([0 1000], [5-sd 5-sd], 'k', 'LineWidth', 1,   'Color', [0.5 0.5 0.5]); hold on; %-grey
         h28 = plot([0 1000], [5+sd 5+sd], 'k', 'LineWidth', 1,   'Color', [0.5 0.5 0.5]); hold on; %-grey
         if package_id == 1
            hy = ylabel('Positive rate (%)', 'Units', 'normalized');
         else
            hy = '';
         end
         if study_id == 1 && package_id == 2
            htitle = title({package; ' '; study_label});
         elseif study_id == 1
            htitle = title({package; ' '; ' '});
         elseif package_id == 2
            htitle = title({' '; study_label});
         end
         if test_type == 1
            xlim([0.8 (5+1)*3+1+0.2]);
         else
            xlim([0.8 (5+1)*5+1+0.2]);
         end
         %-upper limit for ylim cannot be equal lower limit
         if max_y > 0
            ylim([0 max_y]);
         else
            ylim([0 2]);
         end
         %-controlling distance between y axis title and y axis
         hy_pos    = get(hy, 'Position');
         hy_pos(1) = -0.11;
         set(hy, 'Position', hy_pos);
         axes(ax2);
         set(gca, 'XTick', 1000); %-1000: just some value beyond x-axis
         set(gca, 'XTickLabel', repmat({' '}, 1, 1), 'FontSize', 3);
         if package_id > 1
            set(gca, 'YTickLabel', []);
         end
         set([hy htitle], 'FontSize', 5);
         %-you have to find the two values manually (middle locations of first and last HRF labels)
         if test_type == 1
            from = 0.395;
            to   = 0.990;
         else
            from = 0.335;
            to   = 1.055;
         end
         for HRF_model_id = 1:HRF_models_length
            %-specifying font size
            if test_type == 1
               n  = 3;
               fs = 4.5;
            else
               n  = 5;
               fs = 3.25;
            end
            text(ax1, from+(to-from)*(HRF_model_id-1)/(n-1), 0.36, HRF_models_labels{package_id}{HRF_model_id}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', fs, 'Units', 'normalized');
         end   
      end
      fig_ref = subplot(nrows, 3, nrows*3);
      plot(1);
      set(gca, 'visible', 'off');
      %-buffer values: first value down -> text to the left; second value down -> text down
      hlegend = legendflex([h1 h2 h3 h4 h5 h16], {sprintf('boxcar 12s off + 12s on') sprintf('boxcar 16s off + 16s on') sprintf('boxcar 20s off + 20s on') sprintf('CamCAN event-related design') sprintf('~U(3,6)s off + 0.1s on') 'True experimental design'}, 'ref', fig_ref, 'anchor', [4 8], 'buffer', [-342 250], 'box', 'off', 'FontSize', 4, 'nrow', 2);
   end
   figname = [paper '_pos_rates_' test];
   print_to_svg_to_pdf(figname, path_manage);
end


%%%%%%%%%%%%%%%%%%%%%%%% 1st LEVEL: SPATIAL DISTRIBUTION OF SIGNIFICANT CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%
%-type 1: t-test
%-type 2: F-test
%-for 'parfor', some problems
for test_type   = [1 2]
   if test_type == 1
      test = 't';
   else
      test = 'F';
   end
   %-adjusted from the 'autocorr' study, as for other sizes/values, I encountered problems: some figures were a bit shifted...
   for study_id      = range_studies
      for package_id = range_packages
         study       = studies_parameters.study{study_id};
         package     = packages{package_id};
         %-design_scenario 1: true design
         %-design_scenario 2: wrong design
         for design_scenario   = 1:2
            if design_scenario == 1
               exper_design = exper_designs{exper_designs_exp_id(study_id)};
            else
               exper_design = 'event2';
            end
            figure();
            if test_type == 2
               HRF_models_length = 5;
            elseif strcmp(package, 'AFNI') || strcmp(package, 'SPM')
               HRF_models_length = 3;
            else
               HRF_models_length = 2;
            end
            for HRF_model_id     = 1:HRF_models_length
               %10>5 a trick, so that htitles are within the figure window and the subplots do not move
               subplot(10, 1, HRF_model_id+1);
               HRF_model         = HRF_models_labels{package_id}{HRF_model_id};
               load(['combined_results/sp_dist_joint_' test '_' study '_' package '_' exper_design '.mat']);
               sp_dist_joint_adjusted = zeros([4*dims_MNI(1) dims_MNI(2)]);
               for slice_id = 1:4
                  sp_dist_joint_adjusted((slice_id-1)*dims_MNI(1)+1:slice_id*dims_MNI(1), :) = sp_dist_joint((HRF_model_id-1)*dims_MNI(1)+1:HRF_model_id*dims_MNI(1), :, slices(slice_id));
               end
               if design_scenario == 1
                  imagesc(transpose(sp_dist_joint_adjusted), [0 clims_studies(study_id)]);
               else
                  imagesc(transpose(sp_dist_joint_adjusted), [0 2]);
               end
               colormap gray; axis on; axis image; c = colorbar('FontSize', 6);
               set(gca, 'xticklabel', [], 'yticklabel', [], 'xtick', [], 'ytick', []);
               title(['\textbf{' HRF_model '}'], 'FontSize', 6, 'interpreter', 'latex');
               if HRF_model_id == 1
                  %-location can be estimated after getting position of 'title'
                  fig_title = [package ' with ' test '-test'];
                  if package_id == 1 || package_id == 3
                     text(195.6, -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', 8.5);
                  else
                     text(194.0, -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', 8.5);
                  end
               end
               %-'Position': [left bottom width height]
               pos        = get(gca, 'Position');
               pos_bar    = get(c, 'Position');
               pos_bar(1) = 0.75;
               %-2nd term changed, because of pos(2) change, 3rd and 4th terms fixed as otherwise the size of bars can vary...
               pos_bar    = [pos_bar(1) pos(2)-0.023+0.1 0.0224 0.0595];
               set(c,   'Position', pos_bar);
               set(gca, 'Position', [0.04 pos(2) 0.66 0.21]);
            end
            set(gcf, 'units', 'inches', 'position', [0 0 4.0 ceil(10/7 * 7.5)]);
            if design_scenario == 1
               figname = [paper '_1st_level_sp_dist_' test '_' study '_' package '_true_design'];
            else
               figname = [paper '_1st_level_sp_dist_' test '_' study '_' package '_wrong_design'];
            end
            print (['figures/' figname], '-dpdf');
            system(['pdfcrop figures/' figname '.pdf figures/' figname '.pdf']);
         end
      end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%% 2nd LEVEL: SPATIAL DISTRIBUTION OF SIGNIFICANT CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%
%-adjusted from the 'autocorr' study, as for other sizes/values, I encountered problems: some figures were a bit shifted...
%-an example how to overlay two imagesc pics:
%-https://uk.mathworks.com/matlabcentral/answers/373817-i-need-to-overlay-a-color-map-over-a-gray-scale-image-with-with-the-colorbar-on-the-side
%-unfortunately, I have had problems for subplot...
%-btw, 'caxis' doesn't need to be used, while 'Tmap' is simply the second pic
%
%-type 1: t-test
%-type 2: F-test
parfor test_type = [1 2]
   if test_type  == 1
      test = 't';
   else
      test = 'F';
   end
   %-for 'CRIC checkerboard', group brain mask was without V1, which is why 2nd dataset is not considered
   for study_id       = [1 3 4 5]
      for package_id  = range_packages
         figure();
         study        = studies_parameters.study{study_id};
         package      = packages{package_id};
         exper_design = exper_designs{exper_designs_exp_id(study_id)};
         if test_type == 2
            HRF_models_length = 5;
         elseif strcmp(package, 'AFNI') || strcmp(package, 'SPM')
            HRF_models_length = 3;
         else
            HRF_models_length = 2;
         end
         for HRF_model_id  = 1:HRF_models_length
            HRF_model      = HRF_models{package_id}{HRF_model_id};
            %10>5 a trick, so that htitles are within the figure window and the subplots do not move
            brain_template = rot90(niftiread('brain_template_MNI.nii'), 2);
            brain_template = brain_template / max(brain_template(:));
            cluster_binary_location = [path_output study '/' package '/exper_design_' exper_design '/HRF_' HRF_model '/group_analysis_' test '_test/cluster_binary.nii'];
            if exist(cluster_binary_location, 'file') == 2
               cluster_binary = rot90(niftiread(cluster_binary_location), 2);
            else
               disp(['no cluster_binary at ' cluster_binary_location]);
               cluster_binary = zeros(dims_MNI);
            end
            sp_dist_joint_brain_template = zeros([4*dims_MNI(1) dims_MNI(2)]);
            sp_dist_joint_cluster_binary = zeros([4*dims_MNI(1) dims_MNI(2)]);
            for slice_id = 1:4
               sp_dist_joint_brain_template((slice_id-1)*dims_MNI(1)+1:slice_id*dims_MNI(1), :) = brain_template(:, :, slices(slice_id));
               sp_dist_joint_cluster_binary((slice_id-1)*dims_MNI(1)+1:slice_id*dims_MNI(1), :) = cluster_binary(:, :, slices(slice_id));
            end
            sp_dist_joint_brain_template(sp_dist_joint_cluster_binary==1) = 1.2;
            subplot(10, 1, HRF_model_id+1);
            imagesc(transpose(sp_dist_joint_brain_template), [0 1.2]);
            color_palette        = gray();
            %-adding yellow for sig clusters
            color_palette(65, :) = [1 1 0];
            colormap(color_palette);
            axis on; axis image;
            set(gca, 'xticklabel', [], 'yticklabel', [], 'xtick', [], 'ytick', []);
            title(['\textbf{' HRF_models_labels{package_id}{HRF_model_id} '}'], 'FontSize', 6, 'interpreter', 'latex');
            if HRF_model_id == 1
               %-location can be estimated after getting position of 'title'
               fig_title = [package ' with ' test '-test'];
               if package_id == 1 || package_id == 3
                  text(195.6, -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', 8.5);
               else
                  text(194.0, -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontweight', 'bold', 'FontSize', 8.5);
               end
            end
            %-'Position': [left bottom width height]
            pos     = get(gca, 'Position');
            set(gca, 'Position', [0.04 pos(2) 0.66 0.21]);
            set(gcf, 'units', 'inches', 'position', [0 0 4.0 ceil(10/7 * 7.5)]);
         end
         figname = [paper '_2nd_level_sp_dist_' test '_' study '_' package];
         print (['figures/' figname], '-dpdf');
         system(['pdfcrop figures/' figname '.pdf figures/' figname '.pdf']);
      end
   end
end
