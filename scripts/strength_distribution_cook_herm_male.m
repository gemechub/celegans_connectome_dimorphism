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

%cook_herm
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_AM.csv'));
cook_herm_elec = readtable(fullfile(data_path, 'cook_herm_elec_AM.csv'));
cook_herm_comb = readtable(fullfile(data_path, 'cook_herm_combined_AM.csv'));

% cook_male
cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_AM.csv'));% cook_male_electrical
cook_male_elec = readtable(fullfile(data_path, 'cook_male_elec_AM.csv'));
cook_male_comb = readtable(fullfile(data_path, 'cook_male_combined_AM.csv'));


%% Extract local features
%cook herm features
cook_herm_chem_lf = local_features_extraction(cook_herm_chem,cook_cell_types);
cook_herm_elec_lf = local_features_extraction(cook_herm_elec,cook_cell_types);
cook_herm_comb_lf = local_features_extraction(cook_herm_comb,cook_cell_types);

%cook male features
cook_male_chem_lf = local_features_extraction(cook_male_chem,cook_cell_types);
cook_male_elec_lf = local_features_extraction(cook_male_elec,cook_cell_types);
cook_male_comb_lf = local_features_extraction(cook_male_comb,cook_cell_types);

%% strenree survival (1- empeirical cumulative distribution function)

%Herm emperical ecdf fit
%chem
cook_herm_chem_stren = cook_herm_chem_lf.strength;
cook_herm_chem_stren_ccdf = 1 - ecdf(cook_herm_chem_stren); %cumulative distribution
cook_herm_chem_stren_unq = [0; unique(cook_herm_chem_stren)]; %added 0 to balance size

%elec
cook_herm_elec_stren = cook_herm_elec_lf.strength;
cook_herm_elec_stren_ccdf = 1 - ecdf(cook_herm_elec_stren); %cumulative distribution
cook_herm_elec_stren_unq = [0; unique(cook_herm_elec_stren)]; %added 0 to balance size

%comb
cook_herm_comb_stren = cook_herm_comb_lf.strength;
cook_herm_comb_stren_ccdf = 1 - ecdf(cook_herm_comb_stren); %cumulative distribution
cook_herm_comb_stren_unq = [min(cook_herm_comb_stren)-1; unique(cook_herm_comb_stren)]; %added 0 to balance size

%Male
%chem
cook_male_chem_stren = cook_male_chem_lf.strength;
cook_male_chem_stren_ccdf = 1 - ecdf(cook_male_chem_stren); %cumulative distribution
cook_male_chem_stren_unq = [0; unique(cook_male_chem_stren)]; %added 0 to balance size

%elec
cook_male_elec_stren = cook_male_elec_lf.strength;
cook_male_elec_stren_ccdf = 1 - ecdf(cook_male_elec_stren); %cumulative distribution
cook_male_elec_stren_unq = [0; unique(cook_male_elec_stren)]; %added 0 to balance size

%comb
cook_male_comb_stren = cook_male_comb_lf.strength;
cook_male_comb_stren_ccdf = 1 - ecdf(cook_male_comb_stren); %cumulative distribution
cook_male_comb_stren_unq = [0; unique(cook_male_comb_stren)]; %added 0 to balance size


%% Lognormal distribution fitting
%Herm
[cook_herm_chem_stren_logn_ccdf, cook_herm_chem_stren_logn_ksh, cook_herm_chem_stren_logn_ksp,...
    cook_herm_chem_stren_logn_ff, cook_herm_chem_stren_logn_fp] = fit_distribution(cook_herm_chem_stren, 'logn');

[cook_herm_elec_stren_logn_ccdf, cook_herm_elec_stren_logn_ksh, cook_herm_elec_stren_logn_ksp,...
    cook_herm_elec_stren_logn_ff, cook_herm_elec_stren_logn_fp] = fit_distribution(cook_herm_elec_stren, 'logn');

[cook_herm_comb_stren_logn_ccdf, cook_herm_comb_stren_logn_ksh, cook_herm_comb_stren_logn_ksp,...
    cook_herm_comb_stren_logn_ff, cook_herm_comb_stren_logn_fp] = fit_distribution(cook_herm_comb_stren, 'logn');

%Male
[cook_male_chem_stren_logn_ccdf, cook_male_chem_stren_logn_ksh, cook_male_chem_stren_logn_ksp,...
    cook_male_chem_stren_logn_ff, cook_male_chem_stren_logn_fp] = fit_distribution(cook_male_chem_stren, 'logn');

[cook_male_elec_stren_logn_ccdf, cook_male_elec_stren_logn_ksh, cook_male_elec_stren_logn_ksp,...
    cook_male_elec_stren_logn_ff, cook_male_elec_stren_logn_fp] = fit_distribution(cook_male_elec_stren, 'logn');

[cook_male_comb_stren_logn_ccdf, cook_male_comb_stren_logn_ksh, cook_male_comb_stren_logn_ksp,...
    cook_male_comb_stren_logn_ff, cook_male_comb_stren_logn_fp] = fit_distribution(cook_male_comb_stren, 'logn');

%% Exponential distribution fitting
%Herm
[cook_herm_chem_stren_exp_ccdf, cook_herm_chem_stren_exp_ksh, cook_herm_chem_stren_exp_ksp,...
    cook_herm_chem_stren_exp_ff, cook_herm_chem_stren_exp_fp] = fit_distribution(cook_herm_chem_stren, 'exp');

[cook_herm_elec_stren_exp_ccdf, cook_herm_elec_stren_exp_ksh, cook_herm_elec_stren_exp_ksp,...
    cook_herm_elec_stren_exp_ff, cook_herm_elec_stren_exp_fp] = fit_distribution(cook_herm_elec_stren, 'exp');

[cook_herm_comb_stren_exp_ccdf, cook_herm_comb_stren_exp_ksh, cook_herm_comb_stren_exp_ksp,...
    cook_herm_comb_stren_exp_ff, cook_herm_comb_stren_exp_fp] = fit_distribution(cook_herm_comb_stren, 'exp');

%Male
[cook_male_chem_stren_exp_ccdf, cook_male_chem_stren_exp_ksh, cook_male_chem_stren_exp_ksp,...
    cook_male_chem_stren_exp_ff, cook_male_chem_stren_exp_fp] = fit_distribution(cook_male_chem_stren, 'exp');

[cook_male_elec_stren_exp_ccdf, cook_male_elec_stren_exp_ksh, cook_male_elec_stren_exp_ksp,...
    cook_male_elec_stren_exp_ff, cook_male_elec_stren_exp_fp] = fit_distribution(cook_male_elec_stren, 'exp');

[cook_male_comb_stren_exp_ccdf, cook_male_comb_stren_exp_ksh, cook_male_comb_stren_exp_ksp,...
    cook_male_comb_stren_exp_ff, cook_male_comb_stren_exp_fp] = fit_distribution(cook_male_comb_stren, 'exp');

%% weibull distribution fitting
%Herm
[cook_herm_chem_stren_wbl_ccdf, cook_herm_chem_stren_wbl_ksh, cook_herm_chem_stren_wbl_ksp,...
    cook_herm_chem_stren_wbl_ff, cook_herm_chem_stren_wbl_fp] = fit_distribution(cook_herm_chem_stren, 'wbl');

[cook_herm_elec_stren_wbl_ccdf, cook_herm_elec_stren_wbl_ksh, cook_herm_elec_stren_wbl_ksp,...
    cook_herm_elec_stren_wbl_ff, cook_herm_elec_stren_wbl_fp] = fit_distribution(cook_herm_elec_stren, 'wbl');

[cook_herm_comb_stren_wbl_ccdf, cook_herm_comb_stren_wbl_ksh, cook_herm_comb_stren_wbl_ksp,...
    cook_herm_comb_stren_wbl_ff, cook_herm_comb_stren_wbl_fp] = fit_distribution(cook_herm_comb_stren, 'wbl');

%Male
[cook_male_chem_stren_wbl_ccdf, cook_male_chem_stren_wbl_ksh, cook_male_chem_stren_wbl_ksp,...
    cook_male_chem_stren_wbl_ff, cook_male_chem_stren_wbl_fp] = fit_distribution(cook_male_chem_stren, 'wbl');

[cook_male_elec_stren_wbl_ccdf, cook_male_elec_stren_wbl_ksh, cook_male_elec_stren_wbl_ksp,...
    cook_male_elec_stren_wbl_ff, cook_male_elec_stren_wbl_fp] = fit_distribution(cook_male_elec_stren, 'wbl');

[cook_male_comb_stren_wbl_ccdf, cook_male_comb_stren_wbl_ksh, cook_male_comb_stren_wbl_ksp,...
    cook_male_comb_stren_wbl_ff, cook_male_comb_stren_wbl_fp] = fit_distribution(cook_male_comb_stren, 'wbl');


%% Power distribution fitting
%Herm
[cook_herm_chem_stren_pow_ccdf, cook_herm_chem_stren_pow_ksh, cook_herm_chem_stren_pow_ksp,...
    cook_herm_chem_stren_pow_ff, cook_herm_chem_stren_pow_fp] = fit_distribution(cook_herm_chem_stren, 'power');

[cook_herm_elec_stren_pow_ccdf, cook_herm_elec_stren_pow_ksh, cook_herm_elec_stren_pow_ksp,...
    cook_herm_elec_stren_pow_ff, cook_herm_elec_stren_pow_fp] = fit_distribution(cook_herm_elec_stren, 'power');

[cook_herm__comb_stren_pow_ccdf, cook_herm_comb_stren_pow_ksh, cook_herm_comb_stren_pow_ksp,...
    cook_herm_comb_stren_pow_ff, cook_herm_comb_stren_pow_fp] = fit_distribution(cook_herm_comb_stren, 'power');

%Male
[cook_male__chem_stren_pow_ccdf, cook_male_chem_stren_pow_ksh, cook_male_chem_stren_pow_ksp,...
    cook_male_chem_stren_pow_ff, cook_male_chem_stren_pow_fp] = fit_distribution(cook_male_chem_stren, 'power');

[cook_male__elec_stren_pow_ccdf, cook_male_elec_stren_pow_ksh, cook_male_elec_stren_pow_ksp,...
    cook_male_elec_stren_pow_ff, cook_male_elec_stren_pow_fp] = fit_distribution(cook_male_elec_stren, 'power');

[cook_male__comb_stren_pow_ccdf, cook_male_comb_stren_pow_ksh, cook_male_comb_stren_pow_ksp,...
    cook_male_comb_stren_pow_ff, cook_male_comb_stren_pow_fp] = fit_distribution(cook_male_comb_stren, 'power');


%%
fprintf('Herm Chem:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t Exponential:%.4f \t Power:%.4f \t \n',...
    cook_herm_chem_stren_logn_ksp,cook_herm_chem_stren_wbl_ksp,cook_herm_chem_stren_exp_ksp,cook_herm_chem_stren_pow_ksp);


fprintf('Male Chem:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t Exponential:%.4f \t Power:%.4f \t \n\n',...
    cook_male_chem_stren_logn_ksp,cook_male_chem_stren_wbl_ksp,cook_male_chem_stren_exp_ksp,cook_male_chem_stren_pow_ksp);


fprintf('Herm Elec:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t Exponential:%.4f \t Power:%.4f \t \n',...
    cook_herm_elec_stren_logn_ksp,cook_herm_elec_stren_wbl_ksp,cook_herm_elec_stren_exp_ksp,cook_herm_elec_stren_pow_ksp);

fprintf('Male Elec:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t Exponential:%.4f \t Power:%.4f \t \n\n',...
    cook_male_elec_stren_logn_ksp,cook_male_elec_stren_wbl_ksp,cook_male_elec_stren_exp_ksp,cook_male_elec_stren_pow_ksp);


fprintf('Herm Comb:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t Exponential:%.4f \t Power:%.4f \t \n',...
    cook_herm_comb_stren_logn_ksp,cook_herm_comb_stren_wbl_ksp,cook_herm_comb_stren_exp_ksp,cook_herm_comb_stren_pow_ksp);

fprintf('Male Comb:\n');
fprintf('Lognormal: %.4f \t Weibull:%.4f \t Exponential:%.4f \t Power:%.4f \t \n\n',...
    cook_male_comb_stren_logn_ksp,cook_male_comb_stren_wbl_ksp,cook_male_comb_stren_exp_ksp,cook_male_comb_stren_pow_ksp);


%% Comb best fitting model (Logn)
fig = figure;
marker_size = 6;
marker_edge_width = 2;
line_width = 2;

loglog(cook_herm_comb_stren_unq,cook_herm_comb_stren_ccdf,'go','MarkerSize',marker_size);
hold on;
loglog(cook_herm_comb_stren_logn_ff, cook_herm_comb_stren_logn_ccdf,'g','LineWidth',line_width);

loglog(cook_male_comb_stren_unq,cook_male_comb_stren_ccdf,'ro','MarkerSize',marker_size);
loglog(cook_male_comb_stren_logn_ff, cook_male_comb_stren_logn_ccdf,'r','LineWidth',line_width);

lgd = legend('Herm empirical','Herm fitted','Male empirical', 'Male fitted');


xlabel('k'); ylabel('Strength survival function'); 
%set(gca,'FontWeight','bold');

%subplot naming
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;
set(gcf, 'Color', 'W');

set(findobj(gca, 'Type', 'line', 'Marker', 'o'), 'LineWidth', 1.5);

% Remove the top and right axes
ax = gca; % Get current axes
ax.Box = 'off';


% Set the size of the figure
width = 500;  % Width in pixels
height = 350; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_dig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "stren_dist_cook_comb.pdf");
export_fig(fig_save_path);



%% Chem best fitting model (Logn)
fig = figure;

loglog(cook_herm_chem_stren_unq,cook_herm_chem_stren_ccdf,'go','MarkerSize',marker_size);
hold on;
loglog(cook_herm_chem_stren_logn_ff, cook_herm_chem_stren_logn_ccdf,'g','LineWidth',line_width);

loglog(cook_male_chem_stren_unq,cook_male_chem_stren_ccdf,'ro','MarkerSize',marker_size);
loglog(cook_male_chem_stren_logn_ff, cook_male_chem_stren_logn_ccdf,'r','LineWidth',line_width);

lgd = legend('Herm empirical','Herm fitted','Male empirical', 'Male fitted');


xlabel('k'); ylabel('Strength survival function'); 
%set(gca,'FontWeight','bold');

%subplot naming
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;
set(gcf, 'Color', 'W');

set(findobj(gca, 'Type', 'line', 'Marker', 'o'), 'LineWidth', 1.5);

% Remove the top and right axes
ax = gca; % Get current axes
ax.Box = 'off';


% Set the size of the figure
width = 500;  % Width in pixels
height = 350; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_dig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "stren_dist_cook_chem.pdf");
export_fig(fig_save_path);


%% Elec best fitting model (Logn)
fig = figure;

loglog(cook_herm_elec_stren_unq,cook_herm_elec_stren_ccdf,'go','MarkerSize',marker_size);
hold on;
loglog(cook_herm_elec_stren_logn_ff, cook_herm_elec_stren_logn_ccdf,'g','LineWidth',line_width);

loglog(cook_male_elec_stren_unq,cook_male_elec_stren_ccdf,'ro','MarkerSize',marker_size);
loglog(cook_male_elec_stren_logn_ff, cook_male_elec_stren_logn_ccdf,'r','LineWidth',line_width);

lgd = legend('Herm empirical','Herm fitted','Male empirical', 'Male fitted');


xlabel('k'); ylabel('Strength survival function'); 
%set(gca,'FontWeight','bold');

%subplot naming
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;
set(gcf, 'Color', 'W');

set(findobj(gca, 'Type', 'line', 'Marker', 'o'), 'LineWidth', 1.5);

% Remove the top and right axes
ax = gca; % Get current axes
ax.Box = 'off';


% Set the size of the figure
width = 500;  % Width in pixels
height = 350; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_dig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "stren_dist_cook_elec.pdf");
export_fig(fig_save_path);




