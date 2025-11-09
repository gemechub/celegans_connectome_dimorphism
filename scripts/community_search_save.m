%change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")


%% Load datasets

data_path = fullfile(fileparts(pwd), "data");
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_AM.csv'));
cook_herm_elec = readtable(fullfile(data_path, 'cook_herm_elec_AM.csv'));
cook_herm_combined = readtable(fullfile(data_path, 'cook_herm_combined_AM.csv'));

cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_AM.csv'));
cook_male_elec = readtable(fullfile(data_path, 'cook_male_elec_AM.csv'));
cook_male_combined = readtable(fullfile(data_path, 'cook_male_combined_AM.csv'));

%Witvliet dataset
wit7 = readtable(fullfile(data_path, 'wit_d7_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_AM.csv'));

%% communities
%herm
[cook_herm_comb_NC,cook_herm_comb_Q] = community_extraction(cook_herm_combined,1000);
cook_herm_comb_NC_path = fullfile(fileparts(pwd), 'results','cook_herm_comb_communities.csv');
writetable(cook_herm_comb_NC, cook_herm_comb_NC_path);

[cook_herm_chem_NC,cook_herm_chem_Q] = community_extraction(cook_herm_chem,1000);
cook_herm_chem_NC_path = fullfile(fileparts(pwd), 'results','cook_herm_chem_communities.csv');
writetable(cook_herm_chem_NC, cook_herm_chem_NC_path);

[cook_herm_elec_NC,cook_herm_elec_Q] = community_extraction(cook_herm_elec,1000);
cook_herm_elec_NC_path = fullfile(fileparts(pwd), 'results','cook_herm_elec_communities.csv');
writetable(cook_herm_elec_NC, cook_herm_elec_NC_path);

%male
[cook_male_comb_NC,cook_male_comb_Q] = community_extraction(cook_male_combined,1000);
cook_male_comb_NC_path = fullfile(fileparts(pwd), 'results','cook_male_comb_communities.csv');
writetable(cook_male_comb_NC, cook_male_comb_NC_path);

[cook_male_chem_NC,cook_male_chem_Q] = community_extraction(cook_male_chem,1000);
cook_male_chem_NC_path = fullfile(fileparts(pwd), 'results','cook_male_chem_communities.csv');
writetable(cook_male_chem_NC, cook_male_chem_NC_path);

[cook_male_elec_NC,cook_male_elec_Q] = community_extraction(cook_male_elec,1000);
cook_male_elec_NC_path = fullfile(fileparts(pwd), 'results','cook_male_elec_communities.csv');
writetable(cook_male_elec_NC, cook_male_elec_NC_path);

%witv
[wit7_NC,wit7_Q] = community_extraction(wit7,1000);
wit7_NC_path = fullfile(fileparts(pwd), 'results','wit7_communities.csv');
writetable(wit7_NC, wit7_NC_path);

[wit8_NC,wit8_Q] = community_extraction(wit8,1000);
wit8_NC_path = fullfile(fileparts(pwd), 'results','wit8_communities.csv');
writetable(wit8_NC, wit8_NC_path);


%%
clc
fprintf('Modularity and number of communities:\n')
fprintf('\t Herm comb: %d, %d \n', cook_herm_comb_Q, length(unique(cook_herm_comb_NC.Community)));
fprintf('\t Herm chem: %d, %d \n', cook_herm_chem_Q, length(unique(cook_herm_chem_NC.Community)));
fprintf('\t Herm elec: %d, %d \n', cook_herm_elec_Q, length(unique(cook_herm_elec_NC.Community)));

fprintf('\n')
fprintf('\t Male comb: %d, %d \n', cook_male_comb_Q, length(unique(cook_male_comb_NC.Community)));
fprintf('\t Male chem: %d, %d \n', cook_male_chem_Q, length(unique(cook_male_chem_NC.Community)));
fprintf('\t Male elec: %d, %d \n', cook_male_elec_Q, length(unique(cook_male_elec_NC.Community)));

fprintf('\n')
fprintf('\t Wit7: %d, %d \n', wit7_Q, length(unique(wit7_NC.Community)));
fprintf('\t Wit8: %d, %d \n', wit8_Q, length(unique(wit8_NC.Community)));

%save the modularity values
mod_table = table(cook_herm_comb_Q,cook_herm_chem_Q,cook_herm_elec_Q, ...
    cook_male_comb_Q, cook_male_chem_Q, cook_male_elec_Q, wit7_Q, wit8_Q, ...
    {'cook_herm_comb','cook_herm_chem', 'cook_herm_elec', ...
    'cook_male_comb', 'cook_male_chem', 'cook_male_elec', 'wit7', 'wit8'});

mod_path = fullfile(fileparts(pwd), 'results','communities_modularity_values.csv');
writetable(mod_table, mod_path);