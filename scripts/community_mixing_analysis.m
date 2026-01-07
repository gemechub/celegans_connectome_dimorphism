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


cook_herm_chem_el = readtable(fullfile(data_path, 'cook_herm_chem_edgelist.csv'));
cook_herm_elec_el = readtable(fullfile(data_path, 'cook_herm_elec_edgelist.csv'));
cook_herm_combined_el = readtable(fullfile(data_path, 'cook_herm_combined_edgelist.csv'));

cook_male_chem_el = readtable(fullfile(data_path, 'cook_male_chem_edgelist.csv'));
cook_male_elec_el = readtable(fullfile(data_path, 'cook_male_elec_edgelist.csv'));
cook_male_combined_el = readtable(fullfile(data_path, 'cook_male_combined_edgelist.csv'));

%Witvliet dataset
wit7_el = readtable(fullfile(data_path, 'wit_d7_edgelist.csv'));
wit8_el = readtable(fullfile(data_path, 'wit_d8_edgelist.csv'));

%% get mixing matrices
[cook_herm_comb_mixing, ~, cook_herm_comb_mixing_rownorm, ~,cook_herm_comb_mixing_el,...
    cook_herm_comb_K] = community_mixing_matrix(cook_herm_comb_NC, cook_herm_combined_el);


mod_path = fullfile(fileparts(pwd), 'results','cook_herm_combined_community_mixing_matrix.csv');
writetable(cook_herm_comb_mixing, mod_path, 'WriteRowNames', true);

mod_path = fullfile(fileparts(pwd), 'results','cook_herm_combined_community_mixing_matrix_el.csv');
writetable(cook_herm_comb_mixing_el, mod_path, 'WriteRowNames', true);


[cook_male_comb_mixing, ~, cook_male_comb_mixing_rownorm, ~, cook_male_comb_mixing_el,...
    cook_male_comb_K] = community_mixing_matrix(cook_male_comb_NC, cook_male_combined_el);


mod_path = fullfile(fileparts(pwd), 'results','cook_male_combined_community_mixing_matrix.csv');
writetable(cook_male_comb_mixing, mod_path, 'WriteRowNames', true);

mod_path = fullfile(fileparts(pwd), 'results','cook_male_combined_community_mixing_matrix_el.csv');
writetable(cook_male_comb_mixing_el, mod_path, 'WriteRowNames', true);


[wit7_mixing, ~, wit7_mixing_rownorm, ~, wit7_mixing_el,...
    wit7_K] = community_mixing_matrix(wit7_NC, wit7_el);


mod_path = fullfile(fileparts(pwd), 'results','wit7_community_mixing_matrix.csv');
writetable(wit7_mixing, mod_path, 'WriteRowNames', true);

mod_path = fullfile(fileparts(pwd), 'results','wit7_community_mixing_matrix_el.csv');
writetable(wit7_mixing_el, mod_path, 'WriteRowNames', true);


[wit8_mixing, ~, wit8_mixing_rownorm, ~, wit8_mixing_el,...
    wit8_K] = community_mixing_matrix(wit8_NC, wit8_el);


mod_path = fullfile(fileparts(pwd), 'results','wit8_community_mixing_matrix.csv');
writetable(wit8_mixing, mod_path, 'WriteRowNames', true);

mod_path = fullfile(fileparts(pwd), 'results','wit8_community_mixing_matrix_el.csv');
writetable(wit8_mixing_el, mod_path, 'WriteRowNames', true);


%% Simple visualization

fig = figure;
subplot(2,2,1);
colormap(sky);
imagesc(table2array(cook_herm_comb_mixing));
axis image;
colorbar;
title('Herm unnormalized')
xticks(1:cook_herm_comb_K); yticks(1:cook_herm_comb_K);
xlabel('Community');
ylabel('Community');

ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


subplot(2,2,2);
colormap(sky);
imagesc(table2array(cook_herm_comb_mixing_rownorm));
axis image;
colorbar;
title("Herm row normalized");
xticks(1:cook_herm_comb_K); yticks(1:cook_herm_comb_K);
xlabel('Community');
ylabel('Community');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

subplot(2,2,3);
colormap(sky);
imagesc(table2array(cook_male_comb_mixing));
axis image;
colorbar;
title('Male unnormalized')
xticks(1:cook_male_comb_K); yticks(1:cook_male_comb_K);
xlabel('Community');
ylabel('Community');

ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


subplot(2,2,4);
colormap(sky);
imagesc(table2array(cook_male_comb_mixing_rownorm));
axis image;
colorbar;
title("Male row normalized");
xticks(1:cook_male_comb_K); yticks(1:cook_male_comb_K);
xlabel('Community');
ylabel('Community');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


sgtitle("Directed sum of weights within and across communities", 'FontSize', 20);


% Set the size of the figure
width = 800;  % Width in pixels
height = 600; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];


%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "community_mixing_heatmaps.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');

