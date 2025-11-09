%change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% load the datasets

data_path = fullfile(fileparts(pwd), "data");
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

%%
clc
fprintf('Modularity and number of communities:\n')
fprintf('\t Herm comb: %d, %d \n', mod_table.cook_herm_comb_Q, length(unique(cook_herm_comb_NC.Community)));
fprintf('\t Herm chem: %d, %d \n', mod_table.cook_herm_chem_Q, length(unique(cook_herm_chem_NC.Community)));
fprintf('\t Herm elec: %d, %d \n', mod_table.cook_herm_elec_Q, length(unique(cook_herm_elec_NC.Community)));

fprintf('\n')
fprintf('\t Male comb: %d, %d \n', mod_table.cook_male_comb_Q, length(unique(cook_male_comb_NC.Community)));
fprintf('\t Male chem: %d, %d \n', mod_table.cook_male_chem_Q, length(unique(cook_male_chem_NC.Community)));
fprintf('\t Male elec: %d, %d \n', mod_table.cook_male_elec_Q, length(unique(cook_male_elec_NC.Community)));

fprintf('\n')
fprintf('\t Wit7: %d, %d \n', mod_table.wit7_Q, length(unique(wit7_NC.Community)));
fprintf('\t Wit8: %d, %d \n', mod_table.wit8_Q, length(unique(wit8_NC.Community)));


%% count of communities a feature: location
%herm
[cook_herm_chem_unique_communities, cook_herm_chem_unique_locations, cook_herm_chem_location_counts] =...
    community_counts(cook_herm_chem_NC, cook_herm_anatomical_info, 'general_location');

[cook_herm_elec_unique_communities, cook_herm_elec_unique_locations, cook_herm_elec_location_counts] =...
    community_counts(cook_herm_elec_NC, cook_herm_anatomical_info, 'general_location');

[cook_herm_comb_unique_communities, cook_herm_comb_unique_locations, cook_herm_comb_location_counts] =...
    community_counts(cook_herm_comb_NC, cook_herm_anatomical_info, 'general_location');

%male
[cook_male_chem_unique_communities, cook_male_chem_unique_locations, cook_male_chem_location_counts] =...
    community_counts(cook_male_chem_NC, cook_male_anatomical_info, 'general_location');

[cook_male_elec_unique_communities, cook_male_elec_unique_locations, cook_male_elec_location_counts] =...
    community_counts(cook_male_elec_NC, cook_male_anatomical_info, 'general_location');

[cook_male_comb_unique_communities, cook_male_comb_unique_locations, cook_male_comb_location_counts] =...
    community_counts(cook_male_comb_NC, cook_male_anatomical_info, 'general_location');

%witv
[wit7_unique_communities, wit7_unique_locations, wit7_location_counts] =...
    community_counts(wit7_NC, wit_anatomical_info, 'general_location');

[wit8_unique_communities, wit8_unique_locations, wit8_location_counts] =...
    community_counts(wit8_NC, wit_anatomical_info, 'general_location');




%% herm-male general location
fig = figure;
subplot(1,2,1);
[~, cook_herm_comb_location_sorted_com_ind] = sort(sum(cook_herm_comb_location_counts,2),'descend');
cook_herm_comb_location_counts_sorted = cook_herm_comb_location_counts(cook_herm_comb_location_sorted_com_ind,:);
heatmap(cook_herm_comb_unique_locations, cook_herm_comb_unique_communities, cook_herm_comb_location_counts_sorted);
xlabel('General anatomical location')
ylabel('Community');
title('Herm');

ax= gca;
ax.FontName = 'Arial';
ax.FontSize = 16;



subplot(1,2,2);
[~, cook_male_comb_location_sorted_com_ind] = sort(sum(cook_male_comb_location_counts,2),'descend');
cook_male_comb_location_counts_sorted = cook_male_comb_location_counts(cook_male_comb_location_sorted_com_ind,:);
heatmap(cook_male_comb_unique_locations, cook_male_comb_unique_communities, cook_male_comb_location_counts_sorted);
xlabel('General anatomical location')
ylabel('Community');
title('Male');

ax= gca;
ax.FontName = 'Arial';
ax.FontSize = 16;


% Set the size of the figure
width = 900;  % Width in pixels
height = 500; % Height in pixels
set(fig, 'Position', [100, 500, width, height]);

%%
fig_save_path = fullfile(fileparts(pwd), "figure_components", "community_general_location_cook_herm_male.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');


%% cook-witv general location
fig = figure;
subplot(1,3,1);
[~, cook_herm_chem_location_sorted_com_ind] = sort(sum(cook_herm_chem_location_counts,2),'descend');
cook_herm_chem_location_counts_sorted = cook_herm_chem_location_counts(cook_herm_chem_location_sorted_com_ind,:);
heatmap(cook_herm_chem_unique_locations, cook_herm_chem_unique_communities, cook_herm_chem_location_counts_sorted);
xlabel('General anatomical location')
ylabel('Community');
title('Herm');


subplot(1,3,2);
[~, wit7_location_sorted_com_ind] = sort(sum(wit7_location_counts,2),'descend');
wit7_location_counts_sorted = wit7_location_counts(wit7_location_sorted_com_ind,:);
heatmap(wit7_unique_locations, wit7_unique_communities,wit7_location_counts_sorted);
xlabel('General anatomical location')
ylabel('Community');
title('Wit7');

subplot(1,3,3);
[~, wit8_location_sorted_com_ind] = sort(sum(wit8_location_counts,2),'descend');
wit8_location_counts_sorted = wit8_location_counts(wit8_location_sorted_com_ind,:);
heatmap(wit8_unique_locations, wit8_unique_communities,wit8_location_counts_sorted);
xlabel('General anatomical location')
ylabel('Community');
title('Wit8');

% Set the size of the figure
width = 1200;  % Width in pixels
height = 500; % Height in pixels
set(fig, 'Position', [100, 500, width, height]);
%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "community_general_location_cook_wit.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');



