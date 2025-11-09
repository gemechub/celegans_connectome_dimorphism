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
%% degree survival (1- empeirical cumulative distribution function)

%emperical ecdf fit
cook_herm_chem_deg = cook_herm_chem_lf.degree;
cook_herm_chem_deg_ccdf = 1 - ecdf(cook_herm_chem_deg); %cumulative distribution
cook_herm_chem_deg_unq = [0; unique(cook_herm_chem_deg)]; %added 0 to balance size

cook_male_chem_deg = cook_male_chem_lf.degree;
cook_male_chem_deg_ccdf = 1 - ecdf(cook_male_chem_deg); %cumulative distribution
cook_male_chem_deg_unq = [0; unique(cook_male_chem_deg)]; %added 0 to balance size

wit7_deg = wit7_lf.degree;
wit7_deg_ccdf = 1 - ecdf(wit7_deg); 
wit7_deg_unq = [0; unique(wit7_deg)];

wit8_deg = wit8_lf.degree;
wit8_deg_ccdf = 1 - ecdf(wit8_deg); 
wit8_deg_unq = [0; unique(wit8_deg)];

%% lognormal distribution fitting
[cook_herm_chem_deg_logn_ccdf, cook_herm_chem_deg_logn_ksh, cook_herm_chem_deg_logn_ksp,...
    cook_herm_chem_deg_logn_ff, cook_herm_chem_deg_logn_fp] = fit_distribution(cook_herm_chem_deg, 'logn');

[cook_male_chem_deg_logn_ccdf, cook_male_chem_deg_logn_ksh, cook_male_chem_deg_logn_ksp,...
    cook_male_chem_deg_logn_ff, cook_male_chem_deg_logn_fp] = fit_distribution(cook_male_chem_deg, 'logn');

[wit7_deg_logn_ccdf, wit7_deg_logn_ksh, wit7_deg_logn_ksp,wit7_deg_logn_ff, ...
    wit7_deg_logn_fp] = fit_distribution(wit7_deg, 'logn');

[wit8_deg_logn_ccdf, wit8_deg_logn_ksh, wit8_deg_logn_ksp,wit8_deg_logn_ff, ...
    wit8_deg_logn_fp] = fit_distribution(wit8_deg, 'logn');

%% weibull distribution fitting

[cook_herm_chem_deg_wbl_ccdf, cook_herm_chem_deg_wbl_ksh, cook_herm_chem_deg_wbl_ksp,...
    cook_herm_chem_deg_wbl_ff, cook_herm_chem_deg_wbl_fp] = fit_distribution(cook_herm_chem_deg, 'wbl');

[cook_male_chem_deg_wbl_ccdf, cook_male_chem_deg_wbl_ksh, cook_male_chem_deg_wbl_ksp,...
    cook_male_chem_deg_wbl_ff, cook_male_chem_deg_wbl_fp] = fit_distribution(cook_male_chem_deg, 'wbl');

[wit7_deg_wbl_ccdf, wit7_deg_wbl_ksh, wit7_deg_wbl_ksp,wit7_deg_wbl_ff, ...
    wit7_deg_wbl_fp] = fit_distribution(wit7_deg, 'wbl');

[wit8_deg_wbl_ccdf, wit8_deg_wbl_ksh, wit8_deg_wbl_ksp,wit8_deg_wbl_ff, ...
    wit8_deg_wbl_fp] = fit_distribution(wit8_deg, 'wbl');

%%
fprintf('Herm Chem:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    cook_herm_chem_deg_logn_ksp,cook_herm_chem_deg_wbl_ksp);


fprintf('Male Chem:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    cook_male_chem_deg_logn_ksp,cook_male_chem_deg_wbl_ksp);

fprintf('Wit7:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    wit7_deg_logn_ksp,wit7_deg_wbl_ksp);

fprintf('Wit8:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t \n',...
    wit8_deg_logn_ksp,wit8_deg_wbl_ksp);


%% Plots of chem herm - logn , chem male and witv both - weibull
fig = figure;
marker_size = 8;
marker_edge_width = 2;
line_width = 2;
loglog(cook_herm_chem_deg_unq,cook_herm_chem_deg_ccdf,'go', 'MarkerSize',marker_size);
hold on;
loglog(cook_herm_chem_deg_logn_ff, cook_herm_chem_deg_logn_ccdf,'g','LineWidth', line_width);

loglog(cook_male_chem_deg_unq,cook_male_chem_deg_ccdf,'ro','MarkerSize',marker_size);
hold on;
loglog(cook_male_chem_deg_logn_ff, cook_male_chem_deg_wbl_ccdf,'r','LineWidth',line_width);

loglog(wit7_deg_unq, wit7_deg_ccdf, 'bo','MarkerSize',marker_size);
loglog(wit7_deg_logn_ff, wit7_deg_wbl_ccdf,'b','LineWidth',line_width);

loglog(wit8_deg_unq, wit8_deg_ccdf,'o','MarkerSize',marker_size, 'MarkerEdgeColor','magenta');
loglog(wit8_deg_logn_ff, wit8_deg_wbl_ccdf,'magenta','LineWidth',line_width);

xlabel('k'); ylabel('Degree survival function, P(>k)'); 
%set(gca,'FontWeight','bold');
%subplot naming
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;
set(gcf, 'Color', 'W');

set(findobj(gca, 'Type', 'line', 'Marker', 'o'), 'LineWidth', 1.5);

lgd = legend('Cook herm empirical', 'Cook herm fitted','Cook male empirical', 'Cook male fitted',  'Wit7 empricall','Wit7 fitted', 'Wit8 emprical', 'Wit8 fitted');

% Remove the top and right axes
ax = gca; % Get current axes
ax.Box = 'off';

% Set the size of the figure
width = 600;  % Width in pixels
height = 400; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

% Set the position of the legend
% [left, bottom, width, height]
lgd.Position = [0.2, 0.2, 0.2, 0.2];

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_dig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "deg_dist_cook_wit.pdf");
export_fig(fig_save_path);


