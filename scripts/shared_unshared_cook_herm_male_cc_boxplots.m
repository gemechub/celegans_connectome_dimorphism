%% Extract local features
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")


data_path = fullfile(fileparts(pwd), "data");

cook_cell_types = readtable(fullfile(data_path, 'cook_etal_SI4_cell_types.csv'));
wit_cell_types = readtable(fullfile(data_path, 'witvliet_etal_ST1_cell_types.csv'));

%cook_herm_chem
cook_herm_comb = readtable(fullfile(data_path, 'cook_herm_combined_AM.csv'));

% cook_male_chem
cook_male_comb = readtable(fullfile(data_path, 'cook_male_combined_AM.csv'));% cook_male_electrical



%% Extract local features
%cook herm features
cook_herm_comb_lf = local_features_extraction(cook_herm_comb,cook_cell_types);

%cook male features
cook_male_comb_lf = local_features_extraction(cook_male_comb,cook_cell_types);


%% get the sex shared and sex specific neurons' clustering coefficient for box plots
herm_cc = cook_herm_comb_lf.CC_bin;
[herm_shared_cc, herm_specific_cc] = ...
    shared_unshared_features(cook_herm_comb_lf, 'CC_bin', 'herm');

male_cc = cook_male_comb_lf.CC_bin;
[male_shared_cc, male_specific_cc] = ...
    shared_unshared_features(cook_male_comb_lf, 'CC_bin', 'male');

%% 

fig = figure;

subplot(1,2,1);
herm_bp = [herm_cc;herm_shared_cc.CC_bin; herm_specific_cc.CC_bin];
herm_gr = [ones(size(herm_cc));2*ones(size(herm_shared_cc)); 3*ones(size(herm_specific_cc))];
b1 = boxplot(herm_bp, herm_gr, 'Labels',{'All neurons', 'Sex-shared neurons', 'Herm-specific neurons'}, 'Symbol', 'or', 'Width', 0.4);
set(findobj(b1,'Tag','Box'),'LineWidth',2)
set(findobj(b1,'Tag','Median'),'LineWidth',1.5, 'Color', 'r')
set(findobj(b1,'Tag','Upper Whisker'),'LineWidth',1)
set(findobj(b1,'Tag','Lower Whisker'),'LineWidth',1)
ylabel('Clustering coefficient'); 
%set(gca,'FontWeight','bold'); 
%set(gca,'FontWeight','bold');
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 13;
ax.LineWidth = 2;
%ttl = title('A', 'FontSize',18,'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1; 

subplot(1,2,2);
male_bp = [male_cc; male_shared_cc.CC_bin; male_specific_cc.CC_bin];
male_gr = [ones(size(male_cc));2*ones(size(male_shared_cc)); 3*ones(size(male_specific_cc))];
b2 = boxplot(male_bp, male_gr, 'Labels',{'All neurons', 'Sex-shared neurons', 'Male-specific neurons'},'Symbol', 'or', 'Width', 0.4);
set(findobj(b2,'Tag','Box'),'LineWidth',2);
set(findobj(b2,'Tag','Median'),'LineWidth',1.5);
set(findobj(b2,'Tag','Upper Whisker'),'LineWidth',1);
set(findobj(b2,'Tag','Lower Whisker'),'LineWidth',1);
ylabel('Clustering coefficient'); 
%set(gca,'FontWeight','bold');
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 13;
ax.LineWidth = 2;

%ttl = title('B', 'FontSize',18,'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;

% Set the size of the figure
width = 700;  % Width in pixels
height = 300; % Height in pixels
set(fig, 'Position', [100, 500, width, height]);

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_fig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "shared_unshared_cc_cook_herm_male.pdf");
export_fig(fig_save_path);




%% statisitcs 

disp('Herm share vs unshared cc:')
disp(ranksum(herm_shared_cc.CC_bin, herm_specific_cc.CC_bin));

disp('Male share vs unshared degree:')
disp(ranksum(male_shared_cc.CC_bin, male_specific_cc.CC_bin));




