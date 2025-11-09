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

%% similarity of the neurons assigned to the communities
[similarity_score_wit7_wit8, same_diff_matrix_wit7, same_diff_matrix_wit8]... 
    = community_similarity(wit7_NC, wit8_NC);

[similarity_score_wit7_cook, same_diff_matrix_wit7, same_diff_matrix_cook_herm_chem]... 
    = community_similarity(wit7_NC, cook_herm_chem_NC);

[similarity_score_wit8_cook, same_diff_matrix_wit8, same_diff_matrix_cook_herm_chem]... 
    = community_similarity(wit8_NC, cook_herm_chem_NC);

[similarity_score_cook_herm_male, same_diff_matrix_cook_herm_chem, same_diff_matrix_cook_male_chem]...
    =  community_similarity(cook_herm_chem_NC, cook_male_chem_NC);

fprintf('Similarity_score:\n')
fprintf('\t Wit7 vs Wit8: %d \n', similarity_score_wit7_wit8);
fprintf('\t wit7 vs Cook: %d \n', similarity_score_wit7_cook);
fprintf('\t wit8 vs Cook: %d \n', similarity_score_wit8_cook);

fprintf('\t Cook herm vs male: %d \n', similarity_score_cook_herm_male);


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


%% Similarity 
figure;
subplot(2,2,1);
imagesc(same_diff_matrix_cook_herm_chem);
title('Cook herm');
ylabel('Neuron');

subplot(2,2,2);
imagesc(same_diff_matrix_cook_male_chem);
title('Cook male');


subplot(2,2,3);
imagesc(same_diff_matrix_wit7);
title('Wit7');
xlabel('Neuron');
ylabel('Neuron');

subplot(2,2,4);
imagesc(same_diff_matrix_wit8);
title('Wit8');
xlabel('Neuron');


%%
[similarity_score_wit7_wit8, same_diff_matrix_wit7, same_diff_matrix_wit8]... 
    = community_similarity(wit7_NC, wit8_NC);

[similarity_score_wit7_cook, same_diff_matrix_wit7, same_diff_matrix_cook_herm_chem]... 
    = community_similarity(wit7_NC, cook_herm_chem_NC);

[similarity_score_wit8_cook, same_diff_matrix_wit8, same_diff_matrix_cook_herm_chem]... 
    = community_similarity(wit8_NC, cook_herm_chem_NC);



fig = figure();

% % Define custom colormap: rows correspond to [R G B]
cmap = [0 0 0;      % 0 → black
        0 1 0;      % 1 → red
        1 0 0];     % 2 → green

% 
% cmap = [0.2 0.2 0.2;
%         0.8 0.2 0.8;
%         0.0 0.8 0.8];

% cmap = [0 0 0;         % black
%         0.9 0.6 0.0;   % orange
%         0.35 0.7 0.9]; % sky blue

colormap(cmap);
% Fix color limits to match the data values
clim([0 2]);


subplot(2,2,1);
wit7_wit8_diff = 2*same_diff_matrix_wit7-same_diff_matrix_wit8;
wit7_wit8_diff(wit7_wit8_diff==-1)=2;
imagesc(wit7_wit8_diff);
title('Witvliet7 vs Witvliet8', 'FontSize',16);
xlabel('Neuron','FontSize',16);
ylabel('Neuron','FontSize',16);


subplot(2,2,2);
cook_herm_wit7_diff = 2*same_diff_matrix_cook_herm_chem-same_diff_matrix_wit7;
cook_herm_wit7_diff(cook_herm_wit7_diff==-1)=2;
imagesc(cook_herm_wit7_diff);
title('Cook herm vs Witvliet7','FontSize',16);
xlabel('Neuron','FontSize',16);
ylabel('Neuron','FontSize',16);



subplot(2,2,3);
cook_herm_wit8_diff = 2*same_diff_matrix_cook_herm_chem-same_diff_matrix_wit8;
cook_herm_wit8_diff(cook_herm_wit8_diff==-1)=2;
%t1 = ~(t1==0);
imagesc(cook_herm_wit8_diff);
title('Cook herm vs Witvliet8','FontSize',16);
xlabel('Neuron','FontSize',16);
ylabel('Neuron','FontSize',16);



[similarity_score_cook_herm_male, same_diff_matrix_cook_herm_chem, same_diff_matrix_cook_male_chem]...
    =  community_similarity(cook_herm_chem_NC, cook_male_chem_NC);

subplot(2,2,4);
cook_herm_cook_male_diff = 2*same_diff_matrix_cook_herm_chem-same_diff_matrix_cook_male_chem;
cook_herm_cook_male_diff(cook_herm_cook_male_diff==-1)=2;
%t1 = ~(t1==0);
imagesc(cook_herm_cook_male_diff);
title('Cook herm vs Cook male', 'FontSize',16);
xlabel('Neuron','FontSize',16);
ylabel('Neuron','FontSize',16);


cb = colorbar('Position', [0.92 0.1 0.02 0.8], 'Ticks', [0 1 2], 'TickLabels', ...
    {'A','B','C'});
%colorbar('southoutside');


% Add colorbar and label ticks
% cb = colorbar('Position', [0.92 0.1 0.02 0.8], 'Ticks', [0 1 2], 'TickLabels', ...
%     {'In different communities in both connectomes', ...
%     'In the same community in both connectomes', ...
%     'In the same community within only one connectome'});
%cb.Label.String = 'Category';


% Set the size of the figure
width = 800;  % Width in pixels
height = 600; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);


%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "community_similarity_comparison.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');


%%
fig = figure;
similarity_data =[ similarity_score_wit7_wit8,similarity_score_wit7_cook,...
    similarity_score_wit8_cook, similarity_score_cook_herm_male];
bar(similarity_data);
xticklabels({'Witvliet7 vs Witvliet8','Cook herm vs Witvliet7',...
    'Cook herm vs Witvliet8','Cook herm vs Cook male'});
ylabel('Similarity', 'FontSize',18)

ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;
%ttl = title('A', 'FontSize',18,'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1; 

% Set the size of the figure
width = 600;  % Width in pixels
height = 500; % Height in pixels
set(fig, 'Position', [100, 500, width, height]);

%%
fig_save_path = fullfile(fileparts(pwd), "figure_components", "community_similarity_barplot.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');


