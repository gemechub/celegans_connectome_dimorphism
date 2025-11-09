%script to extract local features from the connectivity matrices

clc;
clear;

%change to the directory that contains the script
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% get the main directory, loop through the files and get the CM with its name

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


%% Extract local features
%cook herm features
cook_herm_chem_lf = local_features_extraction(cook_herm_chem,cook_cell_types);
cook_herm_elec_lf = local_features_extraction(cook_herm_elec,cook_cell_types);
cook_herm_comb_lf = local_features_extraction(cook_herm_comb,cook_cell_types);

%cook male features
cook_male_chem_lf = local_features_extraction(cook_male_chem,cook_cell_types);
cook_male_elec_lf = local_features_extraction(cook_male_elec,cook_cell_types);
cook_male_comb_lf = local_features_extraction(cook_male_comb,cook_cell_types);

%witvliet dataset features
wit7_lf = local_features_extraction(wit7,wit_cell_types, 'cell_class');
wit8_lf = local_features_extraction(wit8,wit_cell_types, 'cell_class');

%% Save the local feature tables
result_path = fullfile(fileparts(pwd), "results");

writetable(cook_herm_chem_lf, fullfile(result_path,'cook_herm_chem_local_features.csv'));
writetable(cook_herm_elec_lf, fullfile(result_path,'cook_herm_elec_local_features.csv'));
writetable(cook_herm_comb_lf, fullfile(result_path,'cook_herm_combined_local_features.csv'));

writetable(cook_male_chem_lf, fullfile(result_path,'cook_male_chem_local_features.csv'));
writetable(cook_male_elec_lf, fullfile(result_path,'cook_male_elec_local_features.csv'));
writetable(cook_male_comb_lf, fullfile(result_path,'cook_male_combined_local_features.csv'));

writetable(wit7_lf, fullfile(result_path,'wit7_local_features.csv'));
writetable(wit8_lf, fullfile(result_path,'wit8_local_features.csv'));


   