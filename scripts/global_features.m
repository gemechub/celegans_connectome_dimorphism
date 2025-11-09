%change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% 

data_path = fullfile(fileparts(pwd), "data");
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_strongly_connected_AM.csv'));
cook_herm_elec = readtable(fullfile(data_path, 'cook_herm_elec_strongly_connected_AM.csv'));
cook_herm_combined = readtable(fullfile(data_path, 'cook_herm_combined_strongly_connected_AM.csv'));

cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_strongly_connected_AM.csv'));
cook_male_elec = readtable(fullfile(data_path, 'cook_male_chem_strongly_connected_AM.csv'));
cook_male_combined = readtable(fullfile(data_path, 'cook_male_chem_strongly_connected_AM.csv'));

%Witvliet dataset
wit7 = readtable(fullfile(data_path, 'wit_d7_strongly_connected_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_strongly_connected_AM.csv'));

%% get global features of the connectomes

[gf_bin_herm_chem, gf_wei_herm_chem] = global_features_extraction(cook_herm_chem);
[gf_bin_herm_elec, gf_wei_herm_elec] = global_features_extraction(cook_herm_elec);
[gf_bin_herm_comb, gf_wei_herm_comb] = global_features_extraction(cook_herm_combined);

[gf_bin_male_chem, gf_wei_male_chem] = global_features_extraction(cook_male_chem);
[gf_bin_male_elec, gf_wei_male_elec] = global_features_extraction(cook_male_elec);
[gf_bin_male_comb, gf_wei_male_comb] = global_features_extraction(cook_male_combined);

[gf_bin_wit7, gf_wei_wit7] = global_features_extraction(wit7);
[gf_bin_wit8, gf_wei_wit8] = global_features_extraction(wit8);

%% Save the global feature tables
result_path = fullfile(fileparts(pwd), "results");

%herm binary
writetable(gf_bin_herm_chem, fullfile(result_path,'cook_herm_chem_binary_global_features.csv'));
writetable(gf_bin_herm_elec, fullfile(result_path,'cook_herm_elec_binary_global_features.csv'));
writetable(gf_bin_herm_comb, fullfile(result_path,'cook_herm_combined_binary_global_features.csv'));
%herm weighted
writetable(gf_wei_herm_chem, fullfile(result_path,'cook_herm_chem_weighted_global_features.csv'));
writetable(gf_wei_herm_elec, fullfile(result_path,'cook_herm_elec_weighted_global_features.csv'));
writetable(gf_wei_herm_comb, fullfile(result_path,'cook_herm_weighted_binary_global_features.csv'));
%male binary
writetable(gf_bin_male_chem, fullfile(result_path,'cook_male_chem_binary_global_features.csv'));
writetable(gf_bin_male_elec, fullfile(result_path,'cook_male_elec_binary_global_features.csv'));
writetable(gf_bin_male_comb, fullfile(result_path,'cook_male_combined_binary_global_features.csv'));
%male weighted
writetable(gf_wei_male_chem, fullfile(result_path,'cook_male_chem_weighted_global_features.csv'));
writetable(gf_wei_male_elec, fullfile(result_path,'cook_male_elec_weighted_global_features.csv'));
writetable(gf_wei_male_comb, fullfile(result_path,'cook_male_weighted_binary_global_features.csv'));
%witv7 binary
writetable(gf_bin_wit7, fullfile(result_path,'wit7_binary_global_features.csv'));
%wit7 weighted
writetable(gf_wei_wit7, fullfile(result_path,'wit7_weigted_global_features.csv'));
%witv8 binary
writetable(gf_bin_wit8, fullfile(result_path,'wit8_binary_global_features.csv'));
%wit7 weighted
writetable(gf_wei_wit8, fullfile(result_path,'wit8_weigted_global_features.csv'));



