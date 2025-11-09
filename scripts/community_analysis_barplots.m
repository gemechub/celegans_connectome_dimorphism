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

%% make count of communities feature: cell_type
%herm
[cook_herm_chem_unique_communities, cook_herm_chem_unique_cell_types, cook_herm_chem_cell_type_counts] =...
    community_counts(cook_herm_chem_NC, cook_herm_anatomical_info, 'cell_type');

[cook_herm_elec_unique_communities, cook_herm_elec_unique_cell_types, cook_herm_elec_cell_type_counts] =...
    community_counts(cook_herm_elec_NC, cook_herm_anatomical_info, 'cell_type');

[cook_herm_comb_unique_communities, cook_herm_comb_unique_cell_types, cook_herm_comb_cell_type_counts] =...
    community_counts(cook_herm_comb_NC, cook_herm_anatomical_info, 'cell_type');

%male
[cook_male_chem_unique_communities, cook_male_chem_unique_cell_types, cook_male_chem_cell_type_counts] =...
    community_counts(cook_male_chem_NC, cook_male_anatomical_info, 'cell_type');

[cook_male_elec_unique_communities, cook_male_elec_unique_cell_types, cook_male_elec_cell_type_counts] =...
    community_counts(cook_male_elec_NC, cook_male_anatomical_info, 'cell_type');

[cook_male_comb_unique_communities, cook_male_comb_unique_cell_types, cook_male_comb_cell_type_counts] =...
    community_counts(cook_male_comb_NC, cook_male_anatomical_info, 'cell_type');

%witv
[wit7_unique_communities, wit7_unique_cell_types, wit7_cell_type_counts] =...
    community_counts(wit7_NC, wit_anatomical_info, 'cell_type');

[wit8_unique_communities, wit8_unique_cell_types, wit8_cell_type_counts] =...
    community_counts(wit8_NC, wit_anatomical_info, 'cell_type');

%% count of communities feature: gender
%herm
[~, cook_herm_chem_unique_gender, cook_herm_chem_gender_counts] =...
    community_counts(cook_herm_chem_NC, cook_herm_anatomical_info, 'Gender');

[~, cook_herm_elec_unique_gender, cook_herm_elec_gender_counts] =...
    community_counts(cook_herm_elec_NC, cook_herm_anatomical_info, 'Gender');

[~, cook_herm_comb_unique_gender, cook_herm_comb_gender_counts] =...
    community_counts(cook_herm_comb_NC, cook_herm_anatomical_info, 'Gender');

%male
[~, cook_male_chem_unique_gender, cook_male_chem_gender_counts] =...
    community_counts(cook_male_chem_NC, cook_male_anatomical_info, 'Gender');

[~, cook_male_elec_unique_gender, cook_male_elec_gender_counts] =...
    community_counts(cook_male_elec_NC, cook_male_anatomical_info, 'Gender');

[~, cook_male_comb_unique_gender, cook_male_comb_gender_counts] =...
    community_counts(cook_male_comb_NC, cook_male_anatomical_info, 'Gender');

%witv
[~, wit7_unique_gender, wit7_gender_counts] =...
    community_counts(wit7_NC, wit_anatomical_info, 'Gender');
[~, wit8_unique_gender, wit8_gender_counts] =...
    community_counts(wit8_NC, wit_anatomical_info, 'Gender');



%% Herm-male bar plots
fig = figure();
subplot(2,2,1);

[~, cook_herm_comb_gender_sorted_com_ind] = sort(sum(cook_herm_comb_gender_counts,2),'descend');
cook_herm_comb_gender_counts_sorted = cook_herm_comb_gender_counts(cook_herm_comb_gender_sorted_com_ind,:);
cook_herm_comb_gender_bar_handle = bar(cook_herm_comb_gender_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(cook_herm_comb_unique_gender)); 
for j = 1:numel(cook_herm_comb_unique_gender)
    cook_herm_comb_gender_bar_handle(j).FaceColor = colors(j, :);
end

xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(cook_herm_comb_unique_communities));
xticklabels(cook_herm_comb_unique_communities);
%legend(cook_herm_comb_unique_gender, 'Location', 'northeast');
legend('Sex-shared', 'Hermaphrodite-specific',  'Location', 'northeast');
%ttl = title('A', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

subplot(2,2,2);

[~, cook_male_comb_gender_sorted_com_ind] = sort(sum(cook_male_comb_gender_counts,2),'descend');
cook_male_comb_gender_counts_sorted = cook_male_comb_gender_counts(cook_male_comb_gender_sorted_com_ind,:);
cook_male_comb_gender_bar_handle = bar(cook_male_comb_gender_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(cook_male_comb_unique_gender)); 
for j = 1:numel(cook_male_comb_unique_gender)
    cook_male_comb_gender_bar_handle(j).FaceColor = colors(j, :);
end

xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(cook_male_comb_unique_communities));
xticklabels(cook_male_comb_unique_communities);
%legend(cook_male_comb_unique_gender, 'Location', 'northeast');
legend('Sex-shared', 'Male-specific',  'Location', 'northeast');
%ttl = title('B', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


subplot(2,2,3);
[~, cook_herm_comb_cell_type_sorted_com_ind] = sort(sum(cook_herm_comb_cell_type_counts,2),'descend');
cook_herm_comb_cell_type_counts_sorted = cook_herm_comb_cell_type_counts(cook_herm_comb_cell_type_sorted_com_ind,:);
cook_herm_comb_cell_type_bar_handle = bar(cook_herm_comb_cell_type_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(cook_herm_comb_unique_cell_types)); 
for j = 1:numel(cook_herm_comb_unique_cell_types)
    cook_herm_comb_cell_type_bar_handle(j).FaceColor = colors(j, :);
end

% Set labels and legend
xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(cook_herm_comb_unique_communities));
xticklabels(cook_herm_comb_unique_communities);
%legend(cook_herm_comb_unique_cell_types, 'Location', 'northeast');
legend({'Sensory neurons', 'Interneurons', 'Motor neurons'});

%ttl = title('C', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


subplot(2,2,4);
[~, cook_male_comb_cell_type_sorted_com_ind] = sort(sum(cook_male_comb_cell_type_counts,2),'descend');
cook_male_comb_cell_type_counts_sorted = cook_male_comb_cell_type_counts(cook_male_comb_cell_type_sorted_com_ind,:);
cook_male_comb_cell_type_bar_handle = bar(cook_male_comb_cell_type_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(cook_male_comb_unique_cell_types)); 
for j = 1:numel(cook_male_comb_unique_cell_types)
    cook_male_comb_cell_type_bar_handle(j).FaceColor = colors(j, :);
end

% Set labels and legend
xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(cook_male_comb_unique_communities));
xticklabels(cook_male_comb_unique_communities);
%legend(cook_male_comb_unique_cell_types, 'Location', 'northeast');
legend({'Sensory neurons', 'Interneurons', 'Motor neurons'});

%ttl = title('D', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

% Set the size of the figure
width = 800;  % Width in pixels
height = 600; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];

%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "community_bar_comb_cook.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');





%% cook-witv bar plots
figure();
subplot(2,3,1);

[~, cook_herm_chem_gender_sorted_com_ind] = sort(sum(cook_herm_chem_gender_counts,2),'descend');
cook_herm_chem_gender_counts_sorted = cook_herm_chem_gender_counts(cook_herm_chem_gender_sorted_com_ind,:);
cook_herm_chem_gender_bar_handle = bar(cook_herm_chem_gender_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(cook_herm_chem_unique_gender)); 
for j = 1:numel(cook_herm_chem_unique_gender)
    cook_herm_chem_gender_bar_handle(j).FaceColor = colors(j, :);
end

xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(cook_herm_chem_unique_communities));
xticklabels(cook_herm_chem_unique_communities);
%legend(cook_herm_comb_unique_gender, 'Location', 'northeast');
legend('Sex-shared', 'Hermaphrodite-specific',  'Location', 'northeast');
ttl = title('A', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Times';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;

subplot(2,3,2);
[~, wit7_gender_sorted_com_ind] = sort(sum(wit7_gender_counts,2),'descend');
wit7_gender_counts_sorted = wit7_gender_counts(wit7_gender_sorted_com_ind,:);
wit7_gender_bar_handle = bar(wit7_gender_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(wit7_unique_gender)); 
for j = 1:numel(wit7_unique_gender)
    wit7_gender_bar_handle(j).FaceColor = colors(j, :);
end

xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(wit7_unique_communities));
xticklabels(wit7_unique_communities);
%legend(cook_male_comb_unique_gender, 'Location', 'northeast');
legend('Sex-shared', 'Hermaphrodite-specific',  'Location', 'northeast');
ttl = title('B', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Times';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;



subplot(2,3,3);
[~, wit8_gender_sorted_com_ind] = sort(sum(wit8_gender_counts,2),'descend');
wit8_gender_counts_sorted = wit8_gender_counts(wit8_gender_sorted_com_ind,:);
wit8_gender_bar_handle = bar(wit8_gender_counts_sorted, 'stacked');

% Define colors for each general_location
colors = lines(numel(wit8_unique_gender)); 
for j = 1:numel(wit8_unique_gender)
    wit8_gender_bar_handle(j).FaceColor = colors(j, :);
end
xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(wit8_unique_communities));
xticklabels(wit8_unique_communities);
%legend(cook_male_comb_unique_gender, 'Location', 'northeast');
legend('Sex-shared', 'Hermaphrodite-specific',  'Location', 'northeast');
ttl = title('B', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Times';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;



subplot(2,3,4);
[~, cook_herm_chem_cell_type_sorted_com_ind] = sort(sum(cook_herm_chem_cell_type_counts,2),'descend');
cook_herm_chem_cell_type_counts_sorted = cook_herm_chem_cell_type_counts(cook_herm_chem_cell_type_sorted_com_ind,:);
cook_herm_chem_cell_type_bar_handle = bar(cook_herm_chem_cell_type_counts_sorted, 'stacked');
% Define colors for each general_location
colors = lines(numel(cook_herm_chem_unique_cell_types)); 
for j = 1:numel(cook_herm_chem_unique_cell_types)
    cook_herm_chem_cell_type_bar_handle(j).FaceColor = colors(j, :);
end
% Set labels and legend
xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(cook_herm_chem_unique_communities));
xticklabels(cook_herm_chem_unique_communities);
%legend(cook_herm_comb_unique_cell_types, 'Location', 'northeast');
legend({'Sensory neurons', 'Interneurons', 'Motor neurons'});
ttl = title('C', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Times';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


subplot(2,3,5);
[~, wit7_cell_type_sorted_com_ind] = sort(sum(wit7_cell_type_counts,2),'descend');
wit7_cell_type_counts_sorted = wit7_cell_type_counts(wit7_cell_type_sorted_com_ind,:);
wit7_cell_type_bar_handle = bar(wit7_cell_type_counts_sorted, 'stacked');
% Define colors for each general_location
colors = lines(numel(wit7_unique_cell_types)); 
for j = 1:numel(wit7_unique_cell_types)
    wit7_cell_type_bar_handle(j).FaceColor = colors(j, :);
end
% Set labels and legend
xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(wit7_unique_communities));
xticklabels(wit7_unique_communities);
%legend(cook_herm_comb_unique_cell_types, 'Location', 'northeast');
legend({'Sensory neuron', 'Interneuron', 'Modulatory neuron','Motor neuron'});
ttl = title('C', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Times';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;


subplot(2,3,6);
[~, wit8_cell_type_sorted_com_ind] = sort(sum(wit8_cell_type_counts,2),'descend');
wit8_cell_type_counts_sorted = wit8_cell_type_counts(wit8_cell_type_sorted_com_ind,:);
wit8_cell_type_bar_handle = bar(wit8_cell_type_counts_sorted, 'stacked');
% Define colors for each general_location
colors = lines(numel(wit8_unique_cell_types)); 
for j = 1:numel(wit8_unique_cell_types)
    wit8_cell_type_bar_handle(j).FaceColor = colors(j, :);
end
% Set labels and legend
xlabel('Community');
ylabel('Number of Neurons');
xticks(1:numel(wit8_unique_communities));
xticklabels(wit8_unique_communities);
%legend(cook_herm_comb_unique_cell_types, 'Location', 'northeast');
legend({'Sensory neuron', 'Interneuron', 'Modulatory neuron','Motor neuron'});
ttl = title('C', 'FontSize',18,'FontWeight', 'bold');  ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Times';
ax.FontWeight = 'Bold';
ax.FontSize = 16;
ax.LineWidth = 2;



