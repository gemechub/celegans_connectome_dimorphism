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


%%
cook_herm_chem_fc = functional_carto_features(cook_herm_chem,cook_herm_chem_NC, cook_herm_anatomical_info);
cook_herm_elec_fc = functional_carto_features(cook_herm_elec,cook_herm_elec_NC, cook_herm_anatomical_info);
cook_herm_comb_fc = functional_carto_features(cook_herm_comb,cook_herm_comb_NC, cook_herm_anatomical_info);

cook_male_chem_fc = functional_carto_features(cook_male_chem,cook_male_chem_NC, cook_male_anatomical_info);
cook_male_elec_fc = functional_carto_features(cook_male_elec,cook_male_elec_NC, cook_male_anatomical_info);
cook_male_comb_fc = functional_carto_features(cook_male_comb,cook_male_comb_NC, cook_male_anatomical_info);

wit7_fc = functional_carto_features(wit7,wit7_NC,wit_anatomical_info);
wit8_fc = functional_carto_features(wit8,wit8_NC,wit_anatomical_info);



%% functional cartography: cook herm vs cook male
fig = figure();
subplot(1,2,1);
herm_modz = cook_herm_comb_fc.module_zscore;
herm_partcoef = cook_herm_comb_fc.participation_coefficient;
herm_celltype = cook_herm_comb_fc.cell_type;
herm_celltype_unq = {'sensory', 'interneuron', 'motor'}; 
for j=1:length(herm_celltype_unq)
    curr_celltype = herm_celltype_unq{j};
    ct_idx = find(strcmp(herm_celltype, curr_celltype));
    ct_modz = herm_modz(ct_idx);
    ct_partcoef = herm_partcoef(ct_idx);
    plot(ct_partcoef,ct_modz, '.', 'MarkerSize',15);
    hold on;
end
%set(gca,'FontWeight','bold');
legend(herm_celltype_unq); xlabel('Participation coefficient, P'); ylabel('Within-module z-score, Z'); 
title('Hermaphrodite','FontWeight','normal');

xlim([0 1]); yline(2.5, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
vert_low = linspace(-2, 2.5);
plot(0.05*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.625*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.8*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
vert_up= linspace(2.5,6);
plot(0.3*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.75*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
%ttl = title('A', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax= gca; 
ax.FontName = 'Arial';
set(gcf, 'Color', 'W');

ax.FontSize = 16;
ax.LineWidth = 1.5;


subplot(1,2,2);
male_modz = cook_male_comb_fc.module_zscore;
male_partcoef = cook_male_comb_fc.participation_coefficient;
male_celltype = cook_male_comb_fc.cell_type;
male_celltype_unq = {'sensory', 'interneuron', 'motor'}; 
for j=1:length(male_celltype_unq)
    curr_celltype = male_celltype_unq{j};
    ct_idx = find(strcmp(male_celltype, curr_celltype));
    ct_modz = male_modz(ct_idx);
    ct_partcoef = male_partcoef(ct_idx);
    plot(ct_partcoef,ct_modz, '.', 'MarkerSize',15);
    hold on;
end
%set(gca,'FontWeight','bold');
legend(male_celltype_unq); xlabel('Participation coefficient, P'); ylabel('Within-module z-score, Z'); 
title('Male', 'FontWeight','normal');
xlim([0 1]); yline(2.5, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
vert_low = linspace(-2, 2.5);
plot(0.05*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.625*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.8*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
vert_up= linspace(2.5,7);
plot(0.3*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.75*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');

ax= gca; 
ax.FontName = 'Arial';
%ttl = title('B', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');

ax.FontSize = 16;
ax.LineWidth = 1.5;


% Set the size of the figure
width = 800;  % Width in pixels
height = 400; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];

%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "func_carto_comb_cook.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');


%% functional cartography: wit vs cook
figure();
subplot(1,3,1);
cook_modz = cook_herm_chem_fc.module_zscore;
cook_partcoef = cook_herm_chem_fc.participation_coefficient;
cook_celltype = cook_herm_chem_fc.cell_type;
cook_celltype_unq = {'sensory', 'interneuron', 'motor'}; 
for j=1:length(cook_celltype_unq)
    curr_celltype = cook_celltype_unq{j};
    ct_idx = find(strcmp(cook_celltype, curr_celltype));
    ct_modz = cook_modz(ct_idx);
    ct_partcoef = cook_partcoef(ct_idx);
    plot(ct_partcoef,ct_modz, '.', 'MarkerSize',15);
    hold on;
end
set(gca,'FontWeight','bold');
legend(cook_celltype_unq); xlabel('Participation coefficient, P'); ylabel('Within-module z-score, Z'); title('Cook');
set(gca,'FontWeight','bold');
xlim([0 1]); yline(2.5, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
vert_low = linspace(-2, 2.5);
plot(0.05*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.625*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.8*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
vert_up= linspace(2.5,5);
plot(0.3*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.75*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
ax= gca; ttl = title('A', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');


subplot(1,3,2);
wit7_modz = wit7_fc.module_zscore;
wit7_partcoef = wit7_fc.participation_coefficient;
wit7_celltype = wit7_fc.cell_type;
wit7_celltype_unq = {'Sensory neuron', 'Interneuron', 'Motor neuron', 'Modulatory neuron'}; 
%wit7_celltype_unq = unique(wit7_celltype);
for j=1:length(wit7_celltype_unq)
    curr_celltype = wit7_celltype_unq{j};
    ct_idx = find(strcmp(wit7_celltype, curr_celltype));
    ct_modz = wit7_modz(ct_idx);
    ct_partcoef = wit7_partcoef(ct_idx);
    plot(ct_partcoef,ct_modz, '.', 'MarkerSize',15);
    hold on;
end
set(gca,'FontWeight','bold');
legend(wit7_celltype_unq); xlabel('Participation coefficient, P'); ylabel('Within-module z-score, Z'); title('Male');
set(gca,'FontWeight','bold');
xlim([0 1]); yline(2.5, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
vert_low = linspace(-2, 2.5);
plot(0.05*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.625*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.8*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
vert_up= linspace(2.5,5);
plot(0.3*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.75*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
ax= gca; ttl = title('B', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');

subplot(1,3,3);
wit8_modz = wit8_fc.module_zscore;
wit8_partcoef = wit8_fc.participation_coefficient;
wit8_celltype = wit8_fc.cell_type;
wit8_celltype_unq = {'Sensory neuron', 'Interneuron', 'Motor neuron', 'Modulatory neuron'}; 
%wit8_celltype_unq = unique(wit8_celltype);
for j=1:length(wit8_celltype_unq)
    curr_celltype = wit8_celltype_unq{j};
    ct_idx = find(strcmp(wit8_celltype, curr_celltype));
    ct_modz = wit8_modz(ct_idx);
    ct_partcoef = wit8_partcoef(ct_idx);
    plot(ct_partcoef,ct_modz, '.', 'MarkerSize',15);
    hold on;
end
set(gca,'FontWeight','bold');
legend(wit8_celltype_unq); xlabel('Participation coefficient, P'); ylabel('Within-module z-score, Z'); title('Male');
set(gca,'FontWeight','bold');
xlim([0 1]); yline(2.5, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
vert_low = linspace(-2, 2.5);
plot(0.05*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.625*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.8*ones(size(vert_low)), vert_low, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
vert_up= linspace(2.5,5);
plot(0.3*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
plot(0.75*ones(size(vert_up)), vert_up, 'k', 'LineWidth', 1.5,'HandleVisibility', 'off');
ax= gca; ttl = title('C', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');

%%
herm_chem_region_count = functional_carto_regions_count(cook_herm_chem_fc.module_zscore, cook_herm_chem_fc.participation_coefficient);
herm_elec_region_count = functional_carto_regions_count(cook_herm_elec_fc.module_zscore, cook_herm_elec_fc.participation_coefficient);
herm_comb_region_count = functional_carto_regions_count(cook_herm_comb_fc.module_zscore, cook_herm_comb_fc.participation_coefficient);

male_chem_region_count = functional_carto_regions_count(cook_male_chem_fc.module_zscore, cook_male_chem_fc.participation_coefficient);
male_elec_region_count = functional_carto_regions_count(cook_male_elec_fc.module_zscore, cook_male_elec_fc.participation_coefficient);
male_comb_region_count = functional_carto_regions_count(cook_male_comb_fc.module_zscore, cook_male_comb_fc.participation_coefficient);

wit7_region_count = functional_carto_regions_count(wit7_fc.module_zscore, wit7_fc.participation_coefficient);




