%change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")



%% load the datasets


data_path = fullfile(fileparts(pwd), "data");

cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_strongly_connected_AM.csv'));
cook_herm_elec = readtable(fullfile(data_path, 'cook_herm_elec_strongly_connected_AM.csv'));
cook_herm_combined = readtable(fullfile(data_path, 'cook_herm_combined_strongly_connected_AM.csv'));

cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_strongly_connected_AM.csv'));
cook_male_elec = readtable(fullfile(data_path, 'cook_male_elec_strongly_connected_AM.csv'));
cook_male_combined = readtable(fullfile(data_path, 'cook_male_combined_strongly_connected_AM.csv'));

%Witvliet dataset
wit7 = readtable(fullfile(data_path, 'wit_d7_strongly_connected_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_strongly_connected_AM.csv'));


%null models: herm
cook_herm_chem_null = load(fullfile(data_path, 'cook_herm_chem_sc_null_model.mat'));
cook_herm_chem_null = cook_herm_chem_null.cook_herm_chem_sc_nm;

cook_herm_elec_null = load(fullfile(data_path, 'cook_herm_elec_sc_null_model.mat'));
cook_herm_elec_null = cook_herm_elec_null.cook_herm_elec_sc_nm;

cook_herm_comb_null = load(fullfile(data_path, 'cook_herm_combined_sc_null_model.mat'));
cook_herm_comb_null = cook_herm_comb_null.cook_herm_combined_sc_nm;

%null models: male
cook_male_chem_null = load(fullfile(data_path, 'cook_male_chem_sc_null_model.mat'));
cook_male_chem_null = cook_male_chem_null.cook_male_chem_sc_nm;

cook_male_elec_null = load(fullfile(data_path, 'cook_male_elec_sc_null_model.mat'));
cook_male_elec_null = cook_male_elec_null.cook_male_elec_sc_nm;

cook_male_comb_null = load(fullfile(data_path, 'cook_male_combined_sc_null_model.mat'));
cook_male_comb_null = cook_male_comb_null.cook_male_combined_sc_nm;

%null models: witv
wit7_null = load(fullfile(data_path, 'wit7_sc_null_model.mat'));
wit7_null = wit7_null.wit7_sc_nm;

wit8_null = load(fullfile(data_path, 'wit8_sc_null_model.mat'));
wit8_null = wit8_null.wit8_sc_nm;


%% Global features: cook_herm chem
disp('Cook herm chem:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(cook_herm_chem_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    herm_chem_nulls_bin_gf(i,:) = gf_bin;
    herm_chem_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
herm_chem_bin_globfeat_path = fullfile(fileparts(pwd), 'results','cook_herm_chem_nulls_bin_global_features.csv');
writetable(herm_chem_nulls_bin_gf, herm_chem_bin_globfeat_path);

herm_chem_wei_globfeat_path = fullfile(fileparts(pwd), 'results','cook_herm_chem_nulls_wei_global_features.csv');
writetable(herm_chem_nulls_wei_gf, herm_chem_wei_globfeat_path);

%% Global features: cook_herm elec
disp('Cook herm elec:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(cook_herm_elec_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    herm_elec_nulls_bin_gf(i,:) = gf_bin;
    herm_elec_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
herm_elec_bin_globfeat_path = fullfile(fileparts(pwd), 'results','cook_herm_elec_nulls_bin_global_features.csv');
writetable(herm_elec_nulls_bin_gf, herm_elec_bin_globfeat_path);

herm_elec_wei_globfeat_path = fullfile(fileparts(pwd), 'results','cook_herm_elec_nulls_wei_global_features.csv');
writetable(herm_elec_nulls_wei_gf, herm_elec_wei_globfeat_path);

%% Global features: cook_herm combined
disp('Cook herm combined:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(cook_herm_comb_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    herm_comb_nulls_bin_gf(i,:) = gf_bin;
    herm_comb_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
herm_comb_bin_globfeat_path = fullfile(fileparts(pwd), 'results','cook_herm_combined_nulls_bin_global_features.csv');
writetable(herm_comb_nulls_bin_gf, herm_comb_bin_globfeat_path);

herm_comb_wei_globfeat_path = fullfile(fileparts(pwd), 'results','cook_herm_combined_nulls_wei_global_features.csv');
writetable(herm_comb_nulls_wei_gf, herm_comb_wei_globfeat_path);

%% Global features: cook_male chem
disp('Cook male chem:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(cook_male_chem_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    male_chem_nulls_bin_gf(i,:) = gf_bin;
    male_chem_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
male_chem_bin_globfeat_path = fullfile(fileparts(pwd), 'results','cook_male_chem_nulls_bin_global_features.csv');
writetable(male_chem_nulls_bin_gf, male_chem_bin_globfeat_path);

male_chem_wei_globfeat_path = fullfile(fileparts(pwd), 'results','cook_male_chem_nulls_wei_global_features.csv');
writetable(male_chem_nulls_wei_gf, male_chem_wei_globfeat_path);

%% Global features: cook_male elec
disp('Cook male elec:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(cook_male_elec_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    male_elec_nulls_bin_gf(i,:) = gf_bin;
    male_elec_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
male_elec_bin_globfeat_path = fullfile(fileparts(pwd), 'results','cook_male_elec_nulls_bin_global_features.csv');
writetable(male_elec_nulls_bin_gf, male_elec_bin_globfeat_path);

male_elec_wei_globfeat_path = fullfile(fileparts(pwd), 'results','cook_male_elec_nulls_wei_global_features.csv');
writetable(male_elec_nulls_wei_gf, male_elec_wei_globfeat_path);

%% Global features: cook_male combined
disp('Cook male combined:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(cook_male_comb_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    male_comb_nulls_bin_gf(i,:) = gf_bin;
    male_comb_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
male_comb_bin_globfeat_path = fullfile(fileparts(pwd), 'results','cook_male_combined_nulls_bin_global_features.csv');
writetable(male_comb_nulls_bin_gf, male_comb_bin_globfeat_path);

male_comb_wei_globfeat_path = fullfile(fileparts(pwd), 'results','cook_male_combined_nulls_wei_global_features.csv');
writetable(male_comb_nulls_wei_gf, male_comb_wei_globfeat_path);


%% Global features: wit7
disp('Wit7:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(wit7_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    wit7_nulls_bin_gf(i,:) = gf_bin;
    wit7_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
wit7_bin_globfeat_path = fullfile(fileparts(pwd), 'results','wit7_nulls_bin_global_features.csv');
writetable(wit7_nulls_bin_gf,wit7_bin_globfeat_path);

wit7_wei_globfeat_path = fullfile(fileparts(pwd), 'results','wit7_nulls_wei_global_features.csv');
writetable(wit7_nulls_wei_gf, wit7_wei_globfeat_path);

%% Global features: wit8
disp('Wit8:');
for i=1:1000
    fprintf('%d\n', i);
    cur_null_am = squeeze(wit8_null(i,:,:));
    [gf_bin, gf_wei] = global_features_extraction(cur_null_am,0);
    
    wit8_nulls_bin_gf(i,:) = gf_bin;
    wit8_nulls_wei_gf(i,:) = gf_wei;
end
    
%write tables to csv file
wit8_bin_globfeat_path = fullfile(fileparts(pwd), 'results','wit8_nulls_bin_global_features.csv');
writetable(wit8_nulls_bin_gf,wit8_bin_globfeat_path);

wit8_wei_globfeat_path = fullfile(fileparts(pwd), 'results','wit8_nulls_wei_global_features.csv');
writetable(wit8_nulls_wei_gf, wit8_wei_globfeat_path);


