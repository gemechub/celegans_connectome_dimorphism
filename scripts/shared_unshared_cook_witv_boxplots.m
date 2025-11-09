clc;
clear;
%change to the directory that contains the script
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")


data_path = fullfile(fileparts(pwd), "data");

cook_cell_types = readtable(fullfile(data_path, 'cook_etal_SI4_cell_types.csv'));
wit_cell_types = readtable(fullfile(data_path, 'witvliet_etal_ST1_cell_types.csv'));

%cook_herm_chem
cook_herm_chem = readtable(fullfile(data_path, 'cook_herm_chem_AM.csv'));

% cook_male_chem
cook_male_chem = readtable(fullfile(data_path, 'cook_male_chem_AM.csv'));% cook_male_electrical

%Witvliet
wit7 = readtable(fullfile(data_path, 'wit_d7_AM.csv'));
wit8 = readtable(fullfile(data_path, 'wit_d8_AM.csv'));


%% Extract local features
%cook herm features
cook_herm_chem_lf = local_features_extraction(cook_herm_chem,cook_cell_types);

%cook male features
cook_male_chem_lf = local_features_extraction(cook_male_chem,cook_cell_types);

%witvliet dataset features
wit7_lf = local_features_extraction(wit7,wit_cell_types, 'cell_class');
wit8_lf = local_features_extraction(wit8,wit_cell_types, 'cell_class');


%% 
[cook_herm_chem_shared_degrees, cook_herm_chem_specific_degrees] = ...
    shared_unshared_features(cook_herm_chem_lf, 'degree', 'herm');

[cook_male_chem_shared_degrees, cook_male_chem_specific_degrees] = ...
    shared_unshared_features(cook_male_chem_lf, 'degree', 'male');

[wit7_shared_degrees, wit7_specific_degrees] = ...
    shared_unshared_features(wit7_lf, 'degree', 'herm');

[wit8_shared_degrees, wit8_specific_degrees] = ...
    shared_unshared_features(wit8_lf, 'degree', 'herm');



%% 
fig = figure;

%subplot(1,2,1);
herm_bp = [cook_herm_chem_shared_degrees;cook_male_chem_shared_degrees; wit7_shared_degrees;wit8_shared_degrees];
herm_gr = [ones(size(cook_herm_chem_shared_degrees)); 2*ones(size(cook_male_chem_shared_degrees));...
    3*ones(size(wit7_shared_degrees)); 4*ones(size(wit8_shared_degrees))];
b1 = boxplot(herm_bp.degree, herm_gr, 'Labels',{'Cook herm', 'Cook male', 'Witvliet7', 'Witvliet8'}, 'Symbol', 'or', 'Width', 0.4);
set(findobj(b1,'Tag','Box'),'LineWidth',2)
set(findobj(b1,'Tag','Median'),'LineWidth',1.5, 'Color', 'r')
set(findobj(b1,'Tag','Upper Whisker'),'LineWidth',1)
set(findobj(b1,'Tag','Lower Whisker'),'LineWidth',1)
ylabel('Degree'); 
%set(gca,'FontWeight','bold'); set(gca,'FontWeight','bold');
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 13;
ax.LineWidth = 2;
%ttl = title('A', 'FontSize',18,'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1; 

% Set the size of the figure
width = 350;  % Width in pixels
height = 250; % Height in pixels
set(fig, 'Position', [100, 500, width, height]);

%% exporting the figures
%exporting publication-quality figures is performed using the package
%export_fig: https://github.com/altmany/export_fig 
fig_save_path = fullfile(fileparts(pwd), "figure_components", "shared_unshared_degrees_cook_witv.pdf");
export_fig(fig_save_path);




%% statisitcs 

shared_neurons_degrees = {cook_herm_chem_shared_degrees.degree;cook_male_chem_shared_degrees.degree;...
    wit7_shared_degrees.degree;wit8_shared_degrees.degree};

disp('Cook herm vs Cook male:')
disp(ranksum(cook_herm_chem_shared_degrees.degree,cook_male_chem_shared_degrees.degree));

disp('Cook herm vs wit7:')
disp(ranksum(cook_herm_chem_shared_degrees.degree,wit7_shared_degrees.degree));

disp('Cook herm vs wit8:')
disp(ranksum(cook_herm_chem_shared_degrees.degree,wit8_shared_degrees.degree));

disp('Cook male vs wit7:')
disp(ranksum(cook_male_chem_shared_degrees.degree,wit7_shared_degrees.degree));

disp('Cook male vs wit8:')
disp(ranksum(cook_male_chem_shared_degrees.degree,wit8_shared_degrees.degree));

%%
%wilcoxon rank-sum test if the sex shared and sex-specific degrees have
% Number of vectors
n = length(shared_neurons_degrees);

% Initialize a matrix to store p-values
p_values = zeros(n);

% Perform pairwise ranksum tests
for i = 1:n
    for j = i+1:n
        p_values(i, j) = ranksum(shared_neurons_degrees{i}, shared_neurons_degrees{j});
    end
end

% Display the matrix of p-values
disp('Pairwise ranksum p-values:');
disp(p_values);

%%

shared_neurons_degrees = {cook_herm_chem_shared_degrees.degree; ...
                          cook_male_chem_shared_degrees.degree; ...
                          wit7_shared_degrees.degree; ...
                          wit8_shared_degrees.degree};

% Flatten into one vector with group labels
all_data = [];
group = [];
for i = 1:numel(shared_neurons_degrees)
    all_data = [all_data; shared_neurons_degrees{i}(:)];
    group   = [group; repmat(i, numel(shared_neurons_degrees{i}), 1)];
end

% Nonparametric test: Kruskal-Wallis
[p, tbl, stats] = kruskalwallis(all_data, group, 'off'); % 'off' suppresses the plot

fprintf('Kruskal-Wallis p = %.4f\n', p);

% If significant, do pairwise rank-sum tests
if p < 0.05
    results = multcompare(stats, 'ctype', 'dunn-sidak'); % multiple comparisons correction
    disp('Pairwise comparisons (columns: group1, group2, lowerCI, meanDiff, upperCI, p-value):');
    disp(results);
end


