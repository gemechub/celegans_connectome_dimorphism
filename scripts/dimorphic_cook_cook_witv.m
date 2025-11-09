%change to the directory that contains the script
clear;
clc;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% 

data_path = fullfile(fileparts(pwd), "data");
results_path = fullfile(fileparts(pwd), "results");
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_AM.csv'));
cook_cell_types = readtable(fullfile(data_path, 'cook_etal_SI4_cell_types.csv'));

%Witvliet dataset
wit7 = readtable(fullfile(data_path, 'wit_d7_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_AM.csv'));
wit_cell_types = readtable(fullfile(data_path, 'witvliet_etal_ST1_cell_types.csv'));

%% overlapping neurons local features
cook_witv_lf_overlap = readtable(fullfile(results_path, 'cook_wit_overlapping_neurons_local_features.csv'));

%% dimorphic neurons cook herm vs male
% degree
deg_outlier_cook_wit7 = dimorphic_neurons(cook_witv_lf_overlap, 'degree_cook', 'degree_wit7');
deg_outlier_cook_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'degree_cook', 'degree_wit8');
deg_outlier_wit7_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'degree_wit7', 'degree_wit8');


% strength
stren_outlier_cook_wit7 = dimorphic_neurons(cook_witv_lf_overlap, 'strength_cook', 'strength_wit7');
stren_outlier_cook_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'strength_cook', 'strength_wit8');
stren_outlier_wit7_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'strength_wit7', 'strength_wit8');


% clustering coefficient
cc_outlier_cook_wit7 = dimorphic_neurons(cook_witv_lf_overlap, 'CC_bin_cook', 'CC_bin_wit7');
cc_outlier_cook_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'CC_bin_cook', 'CC_bin_wit8');
cc_outlier_wit7_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'CC_bin_wit7', 'CC_bin_wit8');


% betweenness coefficient
bc_outlier_cook_wit7 = dimorphic_neurons(cook_witv_lf_overlap, 'BC_bin_norm_cook', 'BC_bin_norm_wit7');
bc_outlier_cook_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'BC_bin_norm_cook', 'BC_bin_norm_wit8');
bc_outlier_wit7_wit8 = dimorphic_neurons(cook_witv_lf_overlap, 'BC_bin_norm_wit7', 'BC_bin_norm_wit8');



