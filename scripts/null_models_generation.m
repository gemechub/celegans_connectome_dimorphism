
%change to the directory that contains the script
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% load the tables, extract the AMs and number of nodes

data_path = fullfile(fileparts(pwd), "data");

%original connecitivities
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_AM.csv'));
cook_herm_chem_AM = cook_herm_chem{:,2:end};
cook_herm_chem_node_num = length(cook_herm_chem_AM);

cook_herm_elec = readtable(fullfile(data_path, 'cook_herm_elec_AM.csv'));
cook_herm_elec_AM = cook_herm_elec{:,2:end};
cook_herm_elec_node_num = length(cook_herm_elec_AM);

cook_herm_combined = readtable(fullfile(data_path, 'cook_herm_combined_AM.csv'));
cook_herm_combined_AM = cook_herm_combined{:,2:end};
cook_herm_combined_node_num = length(cook_herm_combined_AM);


cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_AM.csv'));
cook_male_chem_AM = cook_male_chem{:,2:end};
cook_male_chem_node_num = length(cook_male_chem_AM);

cook_male_elec = readtable(fullfile(data_path, 'cook_male_elec_AM.csv'));
cook_male_elec_AM = cook_male_elec{:,2:end};
cook_male_elec_node_num = length(cook_male_elec_AM);

cook_male_combined = readtable(fullfile(data_path, 'cook_male_combined_AM.csv'));
cook_male_combined_AM = cook_male_combined{:,2:end};
cook_male_combined_node_num = length(cook_male_combined_AM);


wit7 = readtable(fullfile(data_path, 'wit_d7_AM.csv'));
wit7_AM = wit7{:,2:end};
wit7_node_num = length(wit7_AM);

wit8 = readtable(fullfile(data_path, 'wit_d8_AM.csv'));
wit8_AM = wit8{:,2:end};
wit8_node_num = length(wit8_AM);

%% fully connected components
%original connecitivities
cook_herm_chem_sc = readtable(fullfile(data_path, 'cook_herm_chem_strongly_connected_AM.csv'));
cook_herm_chem_sc_AM = cook_herm_chem_sc{:,2:end};
cook_herm_chem_sc_node_num = length(cook_herm_chem_sc_AM);

cook_herm_elec_sc = readtable(fullfile(data_path, 'cook_herm_elec_strongly_connected_AM.csv'));
cook_herm_elec_sc_AM = cook_herm_elec_sc{:,2:end};
cook_herm_elec_sc_node_num = length(cook_herm_elec_sc_AM);

cook_herm_combined_sc = readtable(fullfile(data_path, 'cook_herm_combined_strongly_connected_AM.csv'));
cook_herm_combined_sc_AM = cook_herm_combined_sc{:,2:end};
cook_herm_combined_sc_node_num = length(cook_herm_combined_sc_AM);


cook_male_chem_sc = readtable(fullfile(data_path, 'cook_male_chem_strongly_connected_AM.csv'));
cook_male_chem_sc_AM = cook_male_chem_sc{:,2:end};
cook_male_chem_sc_node_num = length(cook_male_chem_sc_AM);

cook_male_elec_sc = readtable(fullfile(data_path, 'cook_male_elec_strongly_connected_AM.csv'));
cook_male_elec_sc_AM = cook_male_elec_sc{:,2:end};
cook_male_elec_sc_node_num = length(cook_male_elec_sc_AM);

cook_male_combined_sc = readtable(fullfile(data_path, 'cook_male_combined_strongly_connected_AM.csv'));
cook_male_combined_sc_AM = cook_male_combined_sc{:,2:end};
cook_male_combined_sc_node_num = length(cook_male_combined_sc_AM);


wit7_sc = readtable(fullfile(data_path, 'wit_d7_strongly_connected_AM.csv'));
wit7_sc_AM = wit7_sc{:,2:end};
wit7_sc_node_num = length(wit7_sc_AM);

wit8_sc = readtable(fullfile(data_path, 'wit_d8_strongly_connected_AM.csv'));
wit8_sc_AM = wit8_sc{:,2:end};
wit8_sc_node_num = length(wit8_sc_AM);


%% null model generation

null_model_num = 1000;

%% cook herm original
cook_herm_chem_nm = zeros(null_model_num, cook_herm_chem_node_num, cook_herm_chem_node_num);
fprintf('Generating model cook_herm_chem:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_herm_chem_nm(i,:,:) = null_model_dir_sign(cook_herm_chem_AM);
end
%save the mat file
save(fullfile(data_path,'cook_herm_chem_null_model.mat'),'cook_herm_chem_nm');

cook_herm_elec_nm = zeros(null_model_num, cook_herm_elec_node_num, cook_herm_elec_node_num);
fprintf('Generating model cook_herm_elec:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_herm_elec_nm(i,:,:) = null_model_dir_sign(cook_herm_elec_AM);
end
%save the mat file
save(fullfile(data_path,'cook_herm_elec_null_model.mat'),'cook_herm_elec_nm');

cook_herm_combined_nm = zeros(null_model_num, cook_herm_combined_node_num, cook_herm_combined_node_num);
fprintf('Generating model cook_herm_combined:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_herm_combined_nm(i,:,:) = null_model_dir_sign(cook_herm_combined_AM);
end
%save the mat file
save(fullfile(data_path,'cook_herm_combined_null_model.mat'),'cook_herm_combined_nm');

%% cook male original

cook_male_chem_nm = zeros(null_model_num, cook_male_chem_node_num, cook_male_chem_node_num);
fprintf('Generating model cook_male_chem:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_male_chem_nm(i,:,:) = null_model_dir_sign(cook_male_chem_AM);
end
%save the mat file
save(fullfile(data_path,'cook_male_chem_null_model.mat'),'cook_male_chem_nm');

cook_male_elec_nm = zeros(null_model_num, cook_male_elec_node_num, cook_male_elec_node_num);
fprintf('Generating model cook_male_elec:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_male_elec_nm(i,:,:) = null_model_dir_sign(cook_male_elec_AM);
end
%save the mat file
save(fullfile(data_path,'cook_male_elec_null_model.mat'),'cook_male_elec_nm');

cook_male_combined_nm = zeros(null_model_num, cook_male_combined_node_num, cook_male_combined_node_num);
fprintf('Generating model cook_male_combined:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_male_combined_nm(i,:,:) = null_model_dir_sign(cook_male_combined_AM);
end
%save the mat file
save(fullfile(data_path,'cook_male_combined_null_model.mat'),'cook_male_combined_nm');

%% Witviliet original

wit7_nm = zeros(null_model_num, wit7_node_num,wit7_node_num);
fprintf('Generating model witv7:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    wit7_nm(i,:,:) = null_model_dir_sign(wit7_AM);
end
%save the mat file
save(fullfile(data_path,'wit7_null_model.mat'),'wit7_nm');


wit8_nm = zeros(null_model_num, wit8_node_num,wit8_node_num);
fprintf('Generating model witv8:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    wit8_nm(i,:,:) = null_model_dir_sign(wit8_AM);

end
%save the mat file
save(fullfile(data_path,'wit8_null_model.mat'),'wit8_nm');

%% cook herm: fully connected

null_model_num = 1000;

cook_herm_chem_sc_nm = zeros(null_model_num, cook_herm_chem_sc_node_num, cook_herm_chem_sc_node_num);
fprintf('Generating model cook_herm_chem_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_herm_chem_sc_nm(i,:,:) = null_model_dir_sign(cook_herm_chem_sc_AM);
end
%save the mat file
save(fullfile(data_path,'cook_herm_chem_sc_null_model.mat'),'cook_herm_chem_sc_nm');

cook_herm_elec_sc_nm = zeros(null_model_num, cook_herm_elec_sc_node_num, cook_herm_elec_sc_node_num);
fprintf('Generating model cook_herm_elec_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_herm_elec_sc_nm(i,:,:) = null_model_dir_sign(cook_herm_elec_sc_AM);
end
%save the mat file
save(fullfile(data_path,'cook_herm_elec_sc_null_model.mat'),'cook_herm_elec_sc_nm');

cook_herm_combined_sc_nm = zeros(null_model_num, cook_herm_combined_sc_node_num, cook_herm_combined_sc_node_num);
fprintf('Generating model cook_herm_combined_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_herm_combined_sc_nm(i,:,:) = null_model_dir_sign(cook_herm_combined_sc_AM);
end
%save the mat file
save(fullfile(data_path,'cook_herm_combined_sc_null_model.mat'),'cook_herm_combined_sc_nm');

%% cook male: fully connected

null_model_num = 1000;

cook_male_chem_sc_nm = zeros(null_model_num, cook_male_chem_sc_node_num, cook_male_chem_sc_node_num);
fprintf('Generating model cook_male_chem_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_male_chem_sc_nm(i,:,:) = null_model_dir_sign(cook_male_chem_sc_AM);
end
%save the mat file
save(fullfile(data_path,'cook_male_chem_sc_null_model.mat'),'cook_male_chem_sc_nm');

cook_male_elec_sc_nm = zeros(null_model_num, cook_male_elec_sc_node_num, cook_male_elec_sc_node_num);
fprintf('Generating model cook_male_elec_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_male_elec_sc_nm(i,:,:) = null_model_dir_sign(cook_male_elec_sc_AM);
end
%save the mat file
save(fullfile(data_path,'cook_male_elec_sc_null_model.mat'),'cook_male_elec_sc_nm');

cook_male_combined_sc_nm = zeros(null_model_num, cook_male_combined_sc_node_num, cook_male_combined_sc_node_num);
fprintf('Generating model cook_male_combined_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    cook_male_combined_sc_nm(i,:,:) = null_model_dir_sign(cook_male_combined_sc_AM);
end
%save the mat file
save(fullfile(data_path,'cook_male_combined_sc_null_model.mat'),'cook_male_combined_sc_nm');


%% Witvliet: fully connected
wit7_sc_nm = zeros(null_model_num, wit7_sc_node_num,wit7_sc_node_num);
fprintf('Generating model witv7_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    wit7_sc_nm(i,:,:) = null_model_dir_sign(wit7_sc_AM);
end
%save the mat file
save(fullfile(data_path,'wit7_sc_null_model.mat'),'wit7_sc_nm');


wit8_sc_nm = zeros(null_model_num, wit8_sc_node_num,wit8_sc_node_num);
fprintf('Generating model witv8_sc:\n');
for i=1:null_model_num
    fprintf('%d\n',i);
    wit8_sc_nm(i,:,:) = null_model_dir_sign(wit8_sc_AM);

end
%save the mat file
save(fullfile(data_path,'wit8_sc_null_model.mat'),'wit8_sc_nm');



