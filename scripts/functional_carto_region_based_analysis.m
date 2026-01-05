%% %change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")


%% load the datasets

%connectivity matrices
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

%Witvliet
wit7 = readtable(fullfile(data_path, 'wit_d7_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_AM.csv'));


%communities
results_path = fullfile(fileparts(pwd), "results");
cook_herm_chem_NC = readtable(fullfile(results_path, 'cook_herm_chem_communities.csv'));
cook_herm_elec_NC = readtable(fullfile(results_path, 'cook_herm_elec_communities.csv'));
cook_herm_comb_NC = readtable(fullfile(results_path, 'cook_herm_comb_communities.csv'));

cook_male_chem_NC = readtable(fullfile(results_path, 'cook_male_chem_communities.csv'));
cook_male_elec_NC = readtable(fullfile(results_path, 'cook_male_elec_communities.csv'));
cook_male_comb_NC = readtable(fullfile(results_path, 'cook_male_comb_communities.csv'));

wit7_NC = readtable(fullfile(results_path, 'wit7_communities.csv'));
wit8_NC = readtable(fullfile(results_path, 'wit8_communities.csv'));

mod_table = readtable(fullfile(results_path, 'communities_modularity_values.csv'));
cook_herm_anatomical_info = readtable(fullfile(data_path, 'cook_herm_neurons_anatomical_info.csv'));
cook_male_anatomical_info = readtable(fullfile(data_path, 'cook_male_neurons_anatomical_info.csv'));
wit_anatomical_info = readtable(fullfile(data_path, 'wit_neurons_anatomical_info.csv'));


%% calculate participation coefficient and within-module z-score
cook_herm_chem_fc = functional_carto_features(cook_herm_chem,cook_herm_chem_NC, cook_herm_anatomical_info);
cook_herm_elec_fc = functional_carto_features(cook_herm_elec,cook_herm_elec_NC, cook_herm_anatomical_info);
cook_herm_comb_fc = functional_carto_features(cook_herm_comb,cook_herm_comb_NC, cook_herm_anatomical_info);

cook_male_chem_fc = functional_carto_features(cook_male_chem,cook_male_chem_NC, cook_male_anatomical_info);
cook_male_elec_fc = functional_carto_features(cook_male_elec,cook_male_elec_NC, cook_male_anatomical_info);
cook_male_comb_fc = functional_carto_features(cook_male_comb,cook_male_comb_NC, cook_male_anatomical_info);

wit7_fc = functional_carto_features(wit7,wit7_NC,wit_anatomical_info);
wit8_fc = functional_carto_features(wit8,wit8_NC,wit_anatomical_info);


%% 
mod_path = fullfile(fileparts(pwd), 'results','cook_herm_chem_functional_cartography_features.csv');
writetable(cook_herm_chem_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','cook_herm_elec_functional_cartography_features.csv');
writetable(cook_herm_elec_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','cook_herm_combined_functional_cartography_features.csv');
writetable(cook_herm_comb_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','cook_male_chem_functional_cartography_features.csv');
writetable(cook_male_chem_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','cook_male_elec_functional_cartography_features.csv');
writetable(cook_male_elec_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','cook_male_combined_functional_cartography_features.csv');
writetable(cook_male_comb_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','wit7_functional_cartography_features.csv');
writetable(wit7_fc, mod_path);

mod_path = fullfile(fileparts(pwd), 'results','wit8_functional_cartography_features.csv');
writetable(wit8_fc, mod_path);


%% perform region based cell-type analysis

% Herm
T_herm_comb = cook_herm_comb_fc;

T_herm_comb.cell_type = categorical(T_herm_comb.cell_type);
T_herm_comb.region    = categorical(T_herm_comb.region);

counts_herm_comb = grpstats(T_herm_comb, {'region','cell_type'}, 'numel', 'DataVars', 'neurons');



region_totals_herm_comb = grpstats(counts_herm_comb, 'region', 'sum', 'DataVars', 'numel_neurons');

counts_herm_comb = outerjoin(counts_herm_comb, ...
    region_totals_herm_comb(:,{'region','sum_numel_neurons'}), ...
    'Keys','region', 'MergeKeys',true);

counts_herm_comb.proportion = counts_herm_comb.numel_neurons ./ counts_herm_comb.sum_numel_neurons;


% Male
T_male_comb = cook_male_comb_fc;

T_male_comb.cell_type = categorical(T_male_comb.cell_type);
T_male_comb.region    = categorical(T_male_comb.region);

counts_male_comb = grpstats(T_male_comb, {'region','cell_type'}, 'numel', 'DataVars', 'neurons');


region_totals_male_comb = grpstats(counts_male_comb, 'region', 'sum', 'DataVars', 'numel_neurons');

counts_male_comb = outerjoin(counts_male_comb, ...
    region_totals_male_comb(:,{'region','sum_numel_neurons'}), ...
    'Keys','region', 'MergeKeys',true);

counts_male_comb.proportion = counts_male_comb.numel_neurons ./ counts_male_comb.sum_numel_neurons;


%% perform region based sex-shared vs sex-specific analysis

% Herm
T_herm_comb = cook_herm_comb_fc;

T_herm_comb.cell_sex = categorical(T_herm_comb.cell_sex);
T_herm_comb.region    = categorical(T_herm_comb.region);

counts_herm_comb_sex = grpstats(T_herm_comb, {'region','cell_sex'}, 'numel', 'DataVars', 'neurons');



region_totals_herm_comb_sex = grpstats(counts_herm_comb_sex, 'region', 'sum', 'DataVars', 'numel_neurons');

counts_herm_comb_sex = outerjoin(counts_herm_comb_sex, ...
    region_totals_herm_comb_sex(:,{'region','sum_numel_neurons'}), ...
    'Keys','region', 'MergeKeys',true);

counts_herm_comb_sex.proportion = counts_herm_comb_sex.numel_neurons ./ counts_herm_comb_sex.sum_numel_neurons;


% Male
T_male_comb = cook_male_comb_fc;

T_male_comb.cell_sex = categorical(T_male_comb.cell_sex);
T_male_comb.region    = categorical(T_male_comb.region);

counts_male_comb_sex = grpstats(T_male_comb, {'region','cell_sex'}, 'numel', 'DataVars', 'neurons');


region_totals_male_comb_sex = grpstats(counts_male_comb_sex, 'region', 'sum', 'DataVars', 'numel_neurons');

counts_male_comb_sex = outerjoin(counts_male_comb_sex, ...
    region_totals_male_comb_sex(:,{'region','sum_numel_neurons'}), ...
    'Keys','region', 'MergeKeys',true);

counts_male_comb_sex.proportion = counts_male_comb_sex.numel_neurons ./ counts_male_comb_sex.sum_numel_neurons;


%% plot the above proportions 

% Proportion tables
prop_table_herm_comb = unstack( ...
    counts_herm_comb(:,{'region','cell_type','proportion'}), ...
    'proportion','cell_type');

prop_table_male_comb = unstack( ...
    counts_male_comb(:,{'region','cell_type','proportion'}), ...
    'proportion','cell_type');

prop_table_herm_comb_sex = unstack( ...
    counts_herm_comb_sex(:,{'region','cell_sex','proportion'}), ...
    'proportion','cell_sex');

prop_table_male_comb_sex = unstack( ...
    counts_male_comb_sex(:,{'region','cell_sex','proportion'}), ...
    'proportion','cell_sex');

% Count tables 
count_table_herm_comb = unstack( ...
    counts_herm_comb(:,{'region','cell_type','numel_neurons'}), ...
    'numel_neurons','cell_type');

count_table_male_comb = unstack( ...
    counts_male_comb(:,{'region','cell_type','numel_neurons'}), ...
    'numel_neurons','cell_type');

count_table_herm_comb_sex = unstack( ...
    counts_herm_comb_sex(:,{'region','cell_sex','numel_neurons'}), ...
    'numel_neurons','cell_sex');

count_table_male_comb_sex = unstack( ...
    counts_male_comb_sex(:,{'region','cell_sex','numel_neurons'}), ...
    'numel_neurons','cell_sex');
%% 
% Proportion
fig = figure();
% Proportions
subplot(2,2,1)
bar(prop_table_herm_comb.region, prop_table_herm_comb{:,2:end}, 'stacked')
ylabel('Proportion')
xlabel('Region')
ylim([0 1])
title('Hermaphrodite (cell type)')
legend(prop_table_herm_comb.Properties.VariableNames(2:end), ...
    'Location','bestoutside', 'Orientation','horizontal')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

subplot(2,2,2)
bar(prop_table_male_comb.region, prop_table_male_comb{:,2:end}, 'stacked')
ylabel('Proportion')
xlabel('Region')
ylim([0 1])
title('Male (cell type)')
legend(prop_table_male_comb.Properties.VariableNames(2:end), ...
     'Location','bestoutside', 'Orientation','horizontal')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off


subplot(2,2,3)
bar(prop_table_herm_comb_sex.region, prop_table_herm_comb_sex{:,2:end}, 'stacked')
ylabel('Proportion')
xlabel('Region')
ylim([0 1])
title('Hermaphrodite (cell sex)')
cellsex_labels  = {'Sex-shared', 'Hermaphrodite-specific'};
legend(cellsex_labels, ...
     'Location','bestoutside', 'Orientation','horizontal')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

subplot(2,2,4)
bar(prop_table_male_comb_sex.region, prop_table_male_comb_sex{:,2:end}, 'stacked')
ylabel('Proportion')
xlabel('Region')
ylim([0 1])
title('Male (cell sex)')

cellsex_labels  = {'Sex-shared', 'Male-specific'};
legend(cellsex_labels, ...
     'Location','bestoutside', 'Orientation','horizontal')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

% figure size
set(fig, 'Position', [100, 100, 1000, 800]);

%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "func_carto_region_cellinfo_proportion.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');



%% Figure : counts
fig = figure();

subplot(2,2,1)
bar(count_table_herm_comb.region, count_table_herm_comb{:,2:end}, 'stacked')
ylabel('Number of neurons')
xlabel('Region')
title('Hermaphrodite (cell type)')
legend(count_table_herm_comb.Properties.VariableNames(2:end), ...
     'Location','northeast', 'Orientation','vertical')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

subplot(2,2,2)
bar(count_table_male_comb.region, count_table_male_comb{:,2:end}, 'stacked')
ylabel('Number of neurons')
xlabel('Region')
title('Male (cell type)')
legend(count_table_male_comb.Properties.VariableNames(2:end), ...
     'Location','northeast', 'Orientation','vertical')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

subplot(2,2,3)
bar(count_table_herm_comb_sex.region, count_table_herm_comb_sex{:,2:end}, 'stacked')
ylabel('Number of neurons')
xlabel('Region')
title('Hermaphrodite (cell sex)')
cellsex_labels  = {'Sex-shared', 'Hermaphrodite-specific'};
legend(cellsex_labels, ...
   'Location','northeast', 'Orientation','vertical')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

subplot(2,2,4)
bar(count_table_male_comb_sex.region, count_table_male_comb_sex{:,2:end}, 'stacked')
ylabel('Number of neurons')
xlabel('Region')
title('Male (cell sex)')
cellsex_labels  = {'Sex-shared', 'Male-specific'};
legend(cellsex_labels, ...
     'Location','northeast', 'Orientation','vertical')
ax = gca; ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=2; box off

set(fig, 'Position', [100, 100, 1000, 700]);

%% 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "func_carto_region_cellinfo_count.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');



%%
% All participation coefficients
pc_herm_all = cook_herm_comb_fc.participation_coefficient;
pc_male_all  = cook_male_comb_fc.participation_coefficient;

% Filtered: sex-shared (GS) neurons only
pc_herm_GS = cook_herm_comb_fc.participation_coefficient( ...
    strcmp(cook_herm_comb_fc.cell_sex,'GS'));
pc_male_GS  = cook_male_comb_fc.participation_coefficient( ...
    strcmp(cook_male_comb_fc.cell_sex,'GS'));

%%
edges = 0:0.05:1;  % fine bins for smooth histogram

fig = figure();

% All neurons subplot: histogram
subplot(1,2,1); hold on
histogram(pc_herm_all, edges, 'Normalization','probability', ...
    'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.5,'EdgeColor','none');
histogram(pc_male_all, edges, 'Normalization','probability', ...
    'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.5,'EdgeColor','none');

xlabel('Participation coefficient')
ylabel('Proportion of neurons')
title('All neurons')
xlim([0 1])
ax = gca; ax.FontName='Arial'; ax.FontSize=14; ax.LineWidth=2; box off

%  GS-only neurons subplot: histogram
subplot(1,2,2); hold on
histogram(pc_herm_GS, edges, 'Normalization','probability', ...
    'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.5,'EdgeColor','none');
histogram(pc_male_GS, edges, 'Normalization','probability', ...
    'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.5,'EdgeColor','none');

xlabel('Participation coefficient')
ylabel('Proportion of neurons')
title('Sex-shared neurons')
xlim([0 1])
ax = gca; ax.FontName='Arial'; ax.FontSize=14; ax.LineWidth=2; box off

% Optional: Adjust figure size
set(fig, 'Position', [100, 100, 1000, 400]);

%% 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "part_coef_histogram.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');

%%
fig = figure();

% All neurons boxplot
subplot(1,2,1)
pc_all = [pc_herm_all; pc_male_all];
group_labels = [repmat({'Hermaphrodite'}, numel(pc_herm_all),1); ...
                repmat({'Male'}, numel(pc_male_all),1)];

h = boxplot(pc_all, group_labels, ...
    'Colors',[0 0.4470 0.7410; 0.8500 0.3250 0.0980], 'Widths',0.5,'Symbol','o');

ax = gca; ax.FontName='Arial'; ax.FontSize=14; ax.LineWidth=2; box off
%ylim([0 1])
ylabel('Participation coefficient')
title('All neurons')

% Make lines thicker
lines = findobj(h,'Type','Line');
for k = 1:length(lines)
    lines(k).LineWidth = 2;
end

hold on
% Overlay points
x_herm = ones(size(pc_herm_all)) + 0.1*(rand(size(pc_herm_all))-0.5);
scatter(x_herm, pc_herm_all, 20, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha',0.5)
x_male = 2*ones(size(pc_male_all)) + 0.1*(rand(size(pc_male_all))-0.5);
scatter(x_male, pc_male_all, 20, [0.8500 0.3250 0.0980], 'filled', 'MarkerFaceAlpha',0.5)

%  GS-only neurons boxplot subplot
subplot(1,2,2)
pc_GS = [pc_herm_GS; pc_male_GS];
group_labels = [repmat({'Hermaphrodite'}, numel(pc_herm_GS),1); ...
                repmat({'Male'}, numel(pc_male_GS),1)];

h = boxplot(pc_GS, group_labels, ...
    'Colors',[0 0.4470 0.7410; 0.8500 0.3250 0.0980], 'Widths',0.5,'Symbol','o');

ax = gca; ax.FontName='Arial'; ax.FontSize=14; ax.LineWidth=2; box off
%ylim([0 1])
ylabel('Participation coefficient')
title('Sex-shared neurons')

% Make lines thicker
lines = findobj(h,'Type','Line');
for k = 1:length(lines)
    lines(k).LineWidth = 2;
end

hold on
x_herm = ones(size(pc_herm_GS)) + 0.1*(rand(size(pc_herm_GS))-0.5);
scatter(x_herm, pc_herm_GS, 20, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha',0.5)
x_male = 2*ones(size(pc_male_GS)) + 0.1*(rand(size(pc_male_GS))-0.5);
scatter(x_male, pc_male_GS, 20, [0.8500 0.3250 0.0980], 'filled', 'MarkerFaceAlpha',0.5)

set(fig, 'Position', [100, 100, 1000, 400]);

%% 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "part_coef_boxplot.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');


%% Statistical tests of the participation coefficients
%Mann-Whitney U test
[p_mw_all,h,stats] = ranksum(pc_herm_all, pc_male_all);
[p_mw_gs,h,stats] = ranksum(pc_herm_GS, pc_male_GS);

[h_ks_all,p_ks_all] = kstest2(pc_herm_all, pc_male_all);
[h_ks_gs,p_ks_gs] = kstest2(pc_herm_GS, pc_male_GS);

%permutation test
nperm = 10000;
diff_obs = median(pc_herm_all) - median(pc_male_all);
all_data = [pc_herm_all; pc_male_all];
n_herm = numel(pc_herm_all);

diff_perm = zeros(nperm,1);
for i = 1:nperm
    perm_idx = randperm(numel(all_data));
    perm_herm = all_data(perm_idx(1:n_herm));
    perm_male = all_data(perm_idx(n_herm+1:end));
    diff_perm(i) = median(perm_herm) - median(perm_male);
end

p_perm_all = mean(abs(diff_perm) >= abs(diff_obs));  % two-sided

%permutation test gs
nperm = 10000;
diff_obs = median(pc_herm_GS) - median(pc_male_GS);
all_data = [pc_herm_GS; pc_male_GS];
n_herm = numel(pc_herm_GS);

diff_perm = zeros(nperm,1);
for i = 1:nperm
    perm_idx = randperm(numel(all_data));
    perm_herm = all_data(perm_idx(1:n_herm));
    perm_male = all_data(perm_idx(n_herm+1:end));
    diff_perm(i) = median(perm_herm) - median(perm_male);
end

p_perm_gs = mean(abs(diff_perm) >= abs(diff_obs));  % two-sided