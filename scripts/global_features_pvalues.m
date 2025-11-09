%change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% load the global features of both original connectomes and their null models
result_path = fullfile(fileparts(pwd), "results");
%herm binary
% original
bin_herm_chem = readtable(fullfile(result_path,'cook_herm_chem_binary_global_features.csv'));
bin_herm_elec = readtable(fullfile(result_path,'cook_herm_elec_binary_global_features.csv'));
bin_herm_comb = readtable(fullfile(result_path,'cook_herm_combined_binary_global_features.csv'));
% null
bin_herm_chem_null = readtable(fullfile(result_path,'cook_herm_chem_nulls_bin_global_features.csv'));
bin_herm_elec_null = readtable(fullfile(result_path,'cook_herm_elec_nulls_bin_global_features.csv'));
bin_herm_comb_null = readtable(fullfile(result_path,'cook_herm_combined_nulls_bin_global_features.csv'));


%herm weighted
% original
wei_herm_chem = readtable(fullfile(result_path,'cook_herm_chem_weighted_global_features.csv'));
wei_herm_elec  = readtable(fullfile(result_path,'cook_herm_elec_weighted_global_features.csv'));
wei_herm_comb = readtable(fullfile(result_path,'cook_herm_combined_weighted_global_features.csv'));
% nulls
wei_herm_chem_null = readtable(fullfile(result_path,'cook_herm_chem_nulls_wei_global_features.csv'));
wei_herm_elec_null  = readtable(fullfile(result_path,'cook_herm_elec_nulls_wei_global_features.csv'));
wei_herm_comb_null = readtable(fullfile(result_path,'cook_herm_combined_nulls_wei_global_features.csv'));


%male binary
% original
bin_male_chem = readtable(fullfile(result_path,'cook_male_chem_binary_global_features.csv'));
bin_male_elec = readtable(fullfile(result_path,'cook_male_elec_binary_global_features.csv'));
bin_male_comb = readtable(fullfile(result_path,'cook_male_combined_binary_global_features.csv'));
% nulls
bin_male_chem_null = readtable(fullfile(result_path,'cook_male_chem_nulls_bin_global_features.csv'));
bin_male_elec_null = readtable(fullfile(result_path,'cook_male_elec_nulls_bin_global_features.csv'));
bin_male_comb_null = readtable(fullfile(result_path,'cook_male_combined_nulls_bin_global_features.csv'));


%male weighted
% original
wei_male_chem = readtable(fullfile(result_path,'cook_male_chem_weighted_global_features.csv'));
wei_male_elec = readtable(fullfile(result_path,'cook_male_elec_weighted_global_features.csv'));
wei_male_comb = readtable(fullfile(result_path,'cook_male_combined_weighted_global_features.csv'));
%nulls
wei_male_chem_null = readtable(fullfile(result_path,'cook_male_chem_nulls_wei_global_features.csv'));
wei_male_elec_null = readtable(fullfile(result_path,'cook_male_elec_nulls_wei_global_features.csv'));
wei_male_comb_null = readtable(fullfile(result_path,'cook_male_combined_nulls_wei_global_features.csv'));

%witv7 binary
% original
bin_wit7 = readtable(fullfile(result_path,'wit7_binary_global_features.csv'));
% nulls
bin_wit7_null = readtable(fullfile(result_path,'wit7_nulls_bin_global_features.csv'));

%wit7 weighted
% original 
wei_wit7 = readtable(fullfile(result_path,'wit7_weigted_global_features.csv'));
% nulls
wei_wit7_null = readtable(fullfile(result_path,'wit7_nulls_wei_global_features.csv'));

%witv8 binary
% original 
bin_wit8 = readtable(fullfile(result_path,'wit8_binary_global_features.csv'));
bin_wit8_null = readtable(fullfile(result_path,'wit8_nulls_bin_global_features.csv'));

%wit7 weighted
% original 
wei_wit8 = readtable(fullfile(result_path,'wit8_weigted_global_features.csv'));
wei_wit8_nulls = readtable(fullfile(result_path,'wit8_nulls_wei_global_features.csv'));

%% statistical test of each connectome with its null model results
%two sided student t-test: bin_herm_chem
[~, cook_herm_chem_pval] = global_feature_ttest(bin_herm_chem_null, bin_herm_chem, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});

[~, cook_herm_elec_pval] = global_feature_ttest(bin_herm_elec_null, bin_herm_elec, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});

[~, cook_herm_comb_pval] = global_feature_ttest(bin_herm_comb_null, bin_herm_comb, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});



[~, cook_male_chem_pval] = global_feature_ttest(bin_male_chem_null, bin_male_chem, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});

[~, cook_male_elec_pval] = global_feature_ttest(bin_male_elec_null, bin_male_elec, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});

[~, cook_male_comb_pval] = global_feature_ttest(bin_male_comb_null, bin_male_comb, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});


[~, wit7_pval] = global_feature_ttest(bin_wit7_null, bin_wit7, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});

[~, wit8_pval] = global_feature_ttest(bin_wit8_null, bin_wit8, ...
    {'cpl', 'assortavity', 'transitivity', 'global_efficiency', 'mean_cc'});

%%

figure;

% Plot CPL Boxplots
subplot(1, 2, 1); % Left subplot for CPL
plot(bin_herm_comb.cpl, 'ro', 'MarkerSize', 8, 'LineWidth', 3); % Herm CPL point
hold on;
b1=boxplot(bin_herm_comb_null.cpl, 'Colors', 'k', 'Symbol', '','Widths', 0.5);
ylim([2.2,3.0]);
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;



subplot(1,2,2);
plot(bin_male_comb.cpl, 'ro', 'MarkerSize', 8, 'LineWidth', 3); % Herm CPL point
hold on;
b1=boxplot(bin_male_comb_null.cpl, 'Colors', 'k', 'Symbol', '','Widths', 0.5);
ylim([2.2,3.0]);
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

%%
fig_save_path = fullfile(fileparts(pwd), "figure_components", "herm_male_comb_cpl.pdf");
export_fig(fig_save_path);

%%

figure;

% Plot global efficiency Boxplots
subplot(1, 2, 1); 
plot(bin_herm_comb.global_efficiency, 'ro', 'MarkerSize', 8, 'LineWidth', 3); 
hold on;
b1=boxplot(bin_herm_comb_null.global_efficiency, 'Colors', 'k', 'Symbol', '','Widths', 0.5);
ylim([0.35,0.5]);
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;



subplot(1,2,2);
plot(bin_male_comb.global_efficiency, 'ro', 'MarkerSize', 8, 'LineWidth', 3); 
hold on;
b1=boxplot(bin_male_comb_null.global_efficiency, 'Colors', 'k', 'Symbol', '','Widths', 0.5);
ylim([0.35,0.5]);
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

%%
fig_save_path = fullfile(fileparts(pwd), "figure_components", "herm_male_comb_global_efficiency.pdf");
export_fig(fig_save_path);

%%

figure;
subplot(1, 2, 1); % Left subplot for CPL
plot(bin_herm_comb.mean_cc, 'ro', 'MarkerSize', 8, 'LineWidth', 3); % Herm CPL point
hold on;
b1=boxplot(bin_herm_comb_null.mean_cc, 'Colors', 'k', 'Symbol', '','Widths',0.5);
ylim([0,0.4]);
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

subplot(1,2,2);
plot(bin_male_comb.mean_cc, 'ro', 'MarkerSize', 8, 'LineWidth', 3); % Herm CPL point
hold on;
b1=boxplot(bin_male_comb_null.mean_cc, 'Colors', 'k', 'Symbol', '','Widths',0.5);
ylim([0,0.4]);
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

%%
fig_save_path = fullfile(fileparts(pwd), "figure_components", "herm_male_comb_mean_cc.pdf");
export_fig(fig_save_path);

%%

figure;
subplot(1, 2, 1); % Left subplot for CPL
plot(bin_herm_comb.assortavity, 'ro', 'MarkerSize', 8, 'LineWidth', 3); % Herm CPL point
hold on;
b1 = boxplot(bin_herm_comb_null.assortavity, 'Colors', 'k', ...
    'Symbol', '', 'Widths', 0.5);
ylim([-0.1,0.2]);
ylabel('Assortativity'); set(gca,'FontWeight','bold');
set(b1,'LineWidth',1);
set(gcf, 'Color', 'W');
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

subplot(1,2,2);
plot(bin_male_comb.assortavity, 'ro', 'MarkerSize', 8, 'LineWidth', 3); % Herm CPL point
hold on;
b1 = boxplot(bin_male_comb_null.assortavity, 'Colors', 'k', ...
    'Symbol', '', 'Widths', 0.5);
ylabel('Assortativity'); set(gca,'FontWeight','bold');
set(b1,'LineWidth',1);
ylim([-0.1,0.2]);
set(gcf, 'Color', 'W');
% set(findobj(b1,'Tag','Box'),'LineWidth',3);
% set(findobj(b1,'Tag','Median'),'LineWidth',3);
% set(findobj(b1,'Tag','Upper Whisker'),'LineWidth',2.5);
% set(findobj(b1,'Tag','Lower Whisker'),'LineWidth',2.5);
ax= gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "herm_male_comb_assort.pdf");
export_fig(fig_save_path);