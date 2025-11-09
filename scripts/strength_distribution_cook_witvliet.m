%script to extract local features from the connectivity matrices
clc;
clear;
%change to the directory that contains the script
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")


data_path = fullfile(fileparts(pwd), "data");

cook_cell_types = readtable(fullfile(data_path, 'cook_etal_SI4_cell_types.csv'));
wit_cell_types = readtable(fullfile(data_path, 'witvliet_etal_ST1_cell_types.csv'));

%cook_herm_chem
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_AM.csv'));

% cook_male_chem
cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_AM.csv'));% cook_male_electrical

%Witvliet
wit7 = readtable(fullfile(data_path, 'wit_d7_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_AM.csv'));


%% Extract local features
%cook herm features
cook_herm_chem_lf = local_features_extraction(cook_herm_chem,cook_cell_types);

%cook male features
cook_male_chem_lf = local_features_extraction(cook_male_chem,cook_cell_types);

%witvliet dataset features
wit7_lf = local_features_extraction(wit7,wit_cell_types, 'cell_class');
wit8_lf = local_features_extraction(wit8,wit_cell_types, 'cell_class');
%% strenree survival (1- empeirical cumulative distribution function)

%emperical ecdf fit
cook_herm_chem_stren = cook_herm_chem_lf.strength;
cook_herm_chem_stren_ccdf = 1 - ecdf(cook_herm_chem_stren); %cumulative distribution
cook_herm_chem_stren_unq = [0; unique(cook_herm_chem_stren)]; %added 0 to balance size

cook_male_chem_stren = cook_male_chem_lf.strength;
cook_male_chem_stren_ccdf = 1 - ecdf(cook_male_chem_stren); %cumulative distribution
cook_male_chem_stren_unq = [0; unique(cook_male_chem_stren)]; %added 0 to balance size

wit7_stren = wit7_lf.strength;
wit7_stren_ccdf = 1 - ecdf(wit7_stren); 
wit7_stren_unq = [0; unique(wit7_stren)];

wit8_stren = wit8_lf.strength;
wit8_stren_ccdf = 1 - ecdf(wit8_stren); 
wit8_stren_unq = [0; unique(wit8_stren)];

%% lognormal distribution fitting
[cook_herm_chem_stren_logn_ccdf, cook_herm_chem_stren_logn_ksh, cook_herm_chem_stren_logn_ksp,...
    cook_herm_chem_stren_logn_ff, cook_herm_chem_stren_logn_fp] = fit_distribution(cook_herm_chem_stren, 'logn');

[cook_male_chem_stren_logn_ccdf, cook_male_chem_stren_logn_ksh, cook_male_chem_stren_logn_ksp,...
    cook_male_chem_stren_logn_ff, cook_male_chem_stren_logn_fp] = fit_distribution(cook_male_chem_stren, 'logn');

[wit7_stren_logn_ccdf, wit7_stren_logn_ksh, wit7_stren_logn_ksp,wit7_stren_logn_ff, ...
    wit7_stren_logn_fp] = fit_distribution(wit7_stren, 'logn');

[wit8_stren_logn_ccdf, wit8_stren_logn_ksh, wit8_stren_logn_ksp,wit8_stren_logn_ff, ...
    wit8_stren_logn_fp] = fit_distribution(wit8_stren, 'logn');

%% weibull distribution fitting

[cook_herm_chem_stren_wbl_ccdf, cook_herm_chem_stren_wbl_ksh, cook_herm_chem_stren_wbl_ksp,...
    cook_herm_chem_stren_wbl_ff, cook_herm_chem_stren_wbl_fp] = fit_distribution(cook_herm_chem_stren, 'wbl');

[cook_male_chem_stren_wbl_ccdf, cook_male_chem_stren_wbl_ksh, cook_male_chem_stren_wbl_ksp,...
    cook_male_chem_stren_wbl_ff, cook_male_chem_stren_wbl_fp] = fit_distribution(cook_male_chem_stren, 'wbl');

[wit7_stren_wbl_ccdf, wit7_stren_wbl_ksh, wit7_stren_wbl_ksp,wit7_stren_wbl_ff, ...
    wit7_stren_wbl_fp] = fit_distribution(wit7_stren, 'wbl');

[wit8_stren_wbl_ccdf, wit8_stren_wbl_ksh, wit8_stren_wbl_ksp,wit8_stren_wbl_ff, ...
    wit8_stren_wbl_fp] = fit_distribution(wit8_stren, 'wbl');

%%
fprintf('Herm Chem:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    cook_herm_chem_stren_logn_ksp,cook_herm_chem_stren_wbl_ksp);


fprintf('Male Chem:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    cook_male_chem_stren_logn_ksp,cook_male_chem_stren_wbl_ksp);

fprintf('Wit7:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    wit7_stren_logn_ksp,wit7_stren_wbl_ksp);

fprintf('Wit8:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    wit8_stren_logn_ksp,wit8_stren_wbl_ksp);

%% Plots of chem herm - logn , chem male and witv both - weibull
fig = figure;
marker_size = 8;
marker_edge_width = 2;
line_width = 2;
loglog(cook_herm_chem_stren_unq,cook_herm_chem_stren_ccdf,'go', 'MarkerSize',marker_size);
hold on;
loglog(cook_herm_chem_stren_logn_ff, cook_herm_chem_stren_logn_ccdf,'g','LineWidth', line_width);

loglog(cook_male_chem_stren_unq,cook_male_chem_stren_ccdf,'ro','MarkerSize',marker_size);
hold on;
loglog(cook_male_chem_stren_logn_ff, cook_male_chem_stren_wbl_ccdf,'r','LineWidth',line_width);

loglog(wit7_stren_unq, wit7_stren_ccdf, 'bo','MarkerSize',marker_size);
loglog(wit7_stren_logn_ff, wit7_stren_wbl_ccdf,'b','LineWidth',line_width);

loglog(wit8_stren_unq, wit8_stren_ccdf,'o','MarkerSize',marker_size, 'MarkerEdgeColor','magenta');
loglog(wit8_stren_logn_ff, wit8_stren_wbl_ccdf,'magenta','LineWidth',line_width);

xlabel('k','FontSize', 22); ylabel('Strength survival function','FontSize', 22); 
%set(gca,'FontWeight','bold');
%subplot naming
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;
set(gcf, 'Color', 'W');

set(findobj(gca, 'Type', 'line', 'Marker', 'o'), 'LineWidth', 1.5);

lgd = legend('Cook herm empirical', 'Cook herm fitted','Cook male empirical', 'Cook male fitted',  'Witvliet7 empricall','Witvliet7 fitted', 'Witvliet8 emprical', 'Witvliet8 fitted');

% Remove the top and right axes
ax = gca; % Get current axes
ax.Box = 'off';

% Set the size of the figure
width = 650;  % Width in pixels
height = 450; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

% Set the position of the legend
% [left, bottom, width, height]
lgd.Position = [0.2, 0.2, 0.2, 0.2];

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_dig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "stren_dist_cook_wit.pdf");
export_fig(fig_save_path);


