%change to the directory that contains the script
clc;
clear;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%add the BCT package to the path
addpath("BCT/")

%% load the datasets


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


%null models: herm
cook_herm_chem_null = load(fullfile(data_path, 'cook_herm_chem_null_model.mat'));
cook_herm_chem_null = cook_herm_chem_null.cook_herm_chem_nm;

cook_herm_elec_null = load(fullfile(data_path, 'cook_herm_elec_null_model.mat'));
cook_herm_elec_null = cook_herm_elec_null.cook_herm_elec_nm;

cook_herm_comb_null = load(fullfile(data_path, 'cook_herm_combined_null_model.mat'));
cook_herm_comb_null = cook_herm_comb_null.cook_herm_combined_nm;

%null models: male
cook_male_chem_null = load(fullfile(data_path, 'cook_male_chem_null_model.mat'));
cook_male_chem_null = cook_male_chem_null.cook_male_chem_nm;

cook_male_elec_null = load(fullfile(data_path, 'cook_male_elec_null_model.mat'));
cook_male_elec_null = cook_male_elec_null.cook_male_elec_nm;

cook_male_comb_null = load(fullfile(data_path, 'cook_male_combined_null_model.mat'));
cook_male_comb_null = cook_male_comb_null.cook_male_combined_nm;

%null models: witv
wit7_null = load(fullfile(data_path, 'wit7_null_model.mat'));
wit7_null = wit7_null.wit7_nm;

wit8_null = load(fullfile(data_path, 'wit8_null_model.mat'));
wit8_null = wit8_null.wit8_nm;

%% local features loading

results_path = fullfile(fileparts(pwd), "results");
cook_herm_chem_lf = readtable(fullfile(results_path, 'cook_herm_chem_local_features.csv'));
cook_herm_elec_lf = readtable(fullfile(results_path, 'cook_herm_elec_local_features.csv'));
cook_herm_combined_lf = readtable(fullfile(results_path, 'cook_herm_combined_local_features.csv'));

cook_male_chem_lf = readtable(fullfile(results_path, 'cook_male_chem_local_features.csv'));
cook_male_elec_lf = readtable(fullfile(results_path, 'cook_male_elec_local_features.csv'));
cook_male_combined_lf = readtable(fullfile(results_path, 'cook_male_combined_local_features.csv'));

%Witvliet dataset
wit7_lf = readtable(fullfile(results_path, 'wit7_local_features.csv'));
wit8_lf = readtable(fullfile(results_path, 'wit8_local_features.csv'));


%% the rich club calculation for original netwokrs and the random networks

%herm chem
disp('Herm chem')
[herm_chem_bin_rcc, herm_chem_bin_rcc_null, herm_chem_norm_bin_rcc, herm_chem_wei_rcc, herm_chem_wei_rcc_null,...
    herm_chem_norm_wei_rcc] = rich_club_extraction(cook_herm_chem,cook_herm_chem_null);

%herm elec
disp('Herm elec')
[herm_elec_bin_rcc, herm_elec_bin_rcc_null, herm_elec_norm_bin_rcc, herm_elec_wei_rcc, herm_elec_wei_rcc_null,...
    herm_elec_norm_wei_rcc] = rich_club_extraction(cook_herm_elec,cook_herm_elec_null);

%herm comb
disp('Herm comb')
[herm_comb_bin_rcc, herm_comb_bin_rcc_null, herm_comb_norm_bin_rcc, herm_comb_wei_rcc, herm_comb_wei_rcc_null,...
    herm_comb_norm_wei_rcc] = rich_club_extraction(cook_herm_combined,cook_herm_comb_null);

%male chem
disp('Male chem')
[male_chem_bin_rcc, male_chem_bin_rcc_null, male_chem_norm_bin_rcc, male_chem_wei_rcc, male_chem_wei_rcc_null,...
    male_chem_norm_wei_rcc] = rich_club_extraction(cook_male_chem,cook_male_chem_null);

%male elec
disp('Male elec')
[male_elec_bin_rcc, male_elec_bin_rcc_null, male_elec_norm_bin_rcc, male_elec_wei_rcc, male_elec_wei_rcc_null,...
    male_elec_norm_wei_rcc] = rich_club_extraction(cook_male_elec,cook_male_elec_null);

%male comb
disp('Male comb')
[male_comb_bin_rcc, male_comb_bin_rcc_null, male_comb_norm_bin_rcc, male_comb_wei_rcc, male_comb_wei_rcc_null,...
    male_comb_norm_wei_rcc] = rich_club_extraction(cook_male_combined,cook_male_comb_null);

%wit7
disp('Wit7')
[wit7_bin_rcc, wit7_bin_rcc_null, wit7_norm_bin_rcc, wit7_wei_rcc, wit7_wei_rcc_null,...
    wit7_norm_wei_rcc] = rich_club_extraction(wit7,wit7_null);

%wit8
disp('Wit8')
[wit8_bin_rcc, wit8_bin_rcc_null, wit8_norm_bin_rcc, wit8_wei_rcc, wit8_wei_rcc_null,...
    wit8_norm_wei_rcc] = rich_club_extraction(wit8,wit8_null);

%% which ones are rich nodes? binary rcc
%herm chem
herm_chem_bin_mean = mean(herm_chem_bin_rcc_null, 'omitnan');
herm_chem_bin_std = std(herm_chem_bin_rcc_null, 'omitnan');
herm_chem_rich_nodes = find((herm_chem_bin_rcc-(herm_chem_bin_mean+3*herm_chem_bin_std))>0.001);
herm_chem_rich_nodes_norm_rcc = herm_chem_norm_bin_rcc(herm_chem_rich_nodes);

%herm elec
herm_elec_bin_mean = mean(herm_elec_bin_rcc_null, 'omitnan');
herm_elec_bin_std = std(herm_elec_bin_rcc_null, 'omitnan');
herm_elec_rich_nodes = find((herm_elec_bin_rcc-(herm_elec_bin_mean+3*herm_elec_bin_std))>0.001);
herm_elec_rich_nodes_norm_rcc = herm_elec_norm_bin_rcc(herm_elec_rich_nodes);

%herm combined
herm_comb_bin_mean = mean(herm_comb_bin_rcc_null, 'omitnan');
herm_comb_bin_std = std(herm_comb_bin_rcc_null, 'omitnan');
herm_comb_rich_nodes = find((herm_comb_bin_rcc-(herm_comb_bin_mean+3*herm_comb_bin_std))>0.001);
herm_comb_rich_nodes_norm_rcc = herm_comb_norm_bin_rcc(herm_comb_rich_nodes);

%male chem
male_chem_bin_mean = mean(male_chem_bin_rcc_null, 'omitnan');
male_chem_bin_std = std(male_chem_bin_rcc_null, 'omitnan');
male_chem_rich_nodes = find((male_chem_bin_rcc-(male_chem_bin_mean+3*male_chem_bin_std))>0.001);
male_chem_rich_nodes_norm_rcc = male_chem_norm_bin_rcc(male_chem_rich_nodes);

%male elec
male_elec_bin_mean = mean(male_elec_bin_rcc_null, 'omitnan');
male_elec_bin_std = std(male_elec_bin_rcc_null, 'omitnan');
male_elec_rich_nodes = find((male_elec_bin_rcc-(male_elec_bin_mean+3*male_elec_bin_std))>0.001);
male_elec_rich_nodes_norm_rcc = male_elec_norm_bin_rcc(male_elec_rich_nodes);

%male combined
male_comb_bin_mean = mean(male_comb_bin_rcc_null, 'omitnan');
male_comb_bin_std = std(male_comb_bin_rcc_null, 'omitnan');
male_comb_rich_nodes = find((male_comb_bin_rcc-(male_comb_bin_mean+3*male_comb_bin_std))>0.001);
male_comb_rich_nodes_norm_rcc = male_comb_norm_bin_rcc(male_comb_rich_nodes);

%wit7
wit7_bin_mean = mean(wit7_bin_rcc_null, 'omitnan');
wit7_bin_std = std(wit7_bin_rcc_null, 'omitnan');
wit7_rich_nodes = find((wit7_bin_rcc-(wit7_bin_mean+3*wit7_bin_std))>0.001);
wit7_rich_nodes_norm_rcc = wit7_norm_bin_rcc(wit7_rich_nodes);

%wit8
wit8_bin_mean = mean(wit8_bin_rcc_null, 'omitnan');
wit8_bin_std = std(wit8_bin_rcc_null, 'omitnan');
wit8_rich_nodes = find((wit8_bin_rcc-(wit8_bin_mean+3*wit8_bin_std))>0.001);
wit8_rich_nodes_norm_rcc = wit8_norm_bin_rcc(wit8_rich_nodes);



%% which ones are rich nodes? weighted rcc

%herm chem
herm_chem_wei_mean = mean(herm_chem_wei_rcc_null, 'omitnan');
herm_chem_wei_std = std(herm_chem_wei_rcc_null, 'omitnan');
herm_chem_wei_rich_nodes = find(herm_chem_norm_wei_rcc>1);
herm_chem_wei_rich_nodes_norm_rcc = herm_chem_norm_wei_rcc(herm_chem_wei_rich_nodes);
herm_chem_wei_rich_neurons = cook_herm_chem_lf.cell_name( ...
    cook_herm_chem_lf.degree>=herm_chem_wei_rich_nodes(1) & ...
    cook_herm_chem_lf.degree<=herm_chem_wei_rich_nodes(end));


%herm elec
herm_elec_wei_mean = mean(herm_elec_wei_rcc_null, 'omitnan');
herm_elec_wei_std = std(herm_elec_wei_rcc_null, 'omitnan');
herm_elec_wei_rich_nodes = find(herm_elec_norm_wei_rcc>1);
herm_elec_wei_rich_nodes_norm_rcc = herm_elec_norm_wei_rcc(herm_elec_wei_rich_nodes);

%herm combined
herm_comb_wei_mean = mean(herm_comb_wei_rcc_null, 'omitnan');
herm_comb_wei_std = std(herm_comb_wei_rcc_null, 'omitnan');
herm_comb_wei_rich_nodes = find(herm_comb_norm_wei_rcc>1);
herm_comb_wei_rich_nodes_norm_rcc = herm_comb_norm_wei_rcc(herm_comb_wei_rich_nodes);


%male chem
male_chem_wei_mean = mean(male_chem_wei_rcc_null, 'omitnan');
male_chem_wei_std = std(male_chem_wei_rcc_null, 'omitnan');
male_chem_wei_rich_nodes = find(male_chem_norm_wei_rcc>1);
male_chem_wei_rich_nodes_norm_rcc = male_chem_norm_wei_rcc(male_chem_wei_rich_nodes);

%male elec
male_elec_wei_mean = mean(male_elec_wei_rcc_null, 'omitnan');
male_elec_wei_std = std(male_elec_wei_rcc_null, 'omitnan');
male_elec_wei_rich_nodes = find(male_elec_norm_wei_rcc>1);
male_elec_wei_rich_nodes_norm_rcc = male_elec_norm_wei_rcc(male_elec_wei_rich_nodes);

%male combined
male_comb_wei_mean = mean(male_comb_wei_rcc_null, 'omitnan');
male_comb_wei_std = std(male_comb_wei_rcc_null, 'omitnan');
male_comb_wei_rich_nodes = find(male_comb_norm_wei_rcc>1);
male_comb_wei_rich_nodes_norm_rcc = male_comb_norm_wei_rcc(male_comb_wei_rich_nodes);

%wit7
wit7_wei_mean = mean(wit7_wei_rcc_null, 'omitnan');
wit7_wei_std = std(wit7_wei_rcc_null, 'omitnan');
wit7_wei_rich_nodes = find(wit7_norm_wei_rcc>1);
wit7_wei_rich_nodes_norm_rcc = wit7_norm_wei_rcc(wit7_wei_rich_nodes);
wit7_wei_rich_neurons = wit7_lf.cell_name( ...
    wit7_lf.degree>=wit7_wei_rich_nodes(1) & wit7_lf.degree<=wit7_wei_rich_nodes(end));

%wit8
wit8_wei_mean = mean(wit8_wei_rcc_null, 'omitnan');
wit8_wei_std = std(wit8_wei_rcc_null, 'omitnan');
wit8_wei_rich_nodes = find(wit8_norm_wei_rcc>1);
wit8_wei_rich_nodes_norm_rcc = wit8_norm_wei_rcc(wit8_wei_rich_nodes);
wit8_wei_rich_neurons = wit8_lf.cell_name( ...
    wit8_lf.degree>=wit8_wei_rich_nodes(1) & wit8_lf.degree<=wit8_wei_rich_nodes(end));

%% comb plots
fig = figure();
%herm unweighted
subplot(2,2,1);
grey  = [220,220,220]./255;
x = [64 93 93 64];
y = [0.02 0.02 2 2];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(herm_comb_bin_rcc,'r','LineWidth',2); 
plot(herm_comb_bin_mean, '-g','LineWidth',2);
plot(herm_comb_norm_bin_rcc,'-b','LineWidth',1.5);
plot(herm_comb_rich_nodes,herm_comb_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K');
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1, 'LineWidth', 1.5);
l = legend('\Phi(k)_{Herm}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
l.Position = [0.375, 0.85, 0.1, 0.1];
legend boxoff;

ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 12;
ax.LineWidth = 2;
ttl = title('A', 'FontSize',14, 'FontName','Arial','FontWeight','normal'); 
ttl.Units = 'Normalize'; ttl.Position(1) = -0.2; ttl.Position(2) = 1.05;
set(gcf, 'Color', 'W');
%title('Unweighted rich club curve(herm)'); 


%male unweighted
subplot(2,2,2);
grey  = [220,220,220]./255;
x = [60 87 87 60];
y = [0.02 0.02 2.5 2.5];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(male_comb_bin_rcc,'-r','LineWidth',2); 
plot(male_comb_bin_mean, '-g','LineWidth',2);
plot(male_comb_norm_bin_rcc,'-b','LineWidth',1.5);
plot(male_comb_rich_nodes,male_comb_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K'); 
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1,'LineWidth', 1.5);
l= legend('\Phi(k)_{Male}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
l.Position = [0.795, 0.85, 0.1, 0.1];
legend boxoff;
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 12;
ax.LineWidth = 2;
ttl = title('B', 'FontSize',14, 'FontName','Arial','FontWeight','normal'); 
ttl.Units = 'Normalize'; ttl.Position(1) = -0.2; ttl.Position(2) = 1.05;

set(gcf, 'Color', 'W');


%herm weighted 
subplot(2,2,3);
plot(herm_comb_wei_rcc,'r','LineWidth',2);
hold on;
plot(herm_comb_wei_mean, '-g','LineWidth',2);
plot(herm_comb_norm_wei_rcc,'b','LineWidth',1.5);
plot(herm_comb_wei_rich_nodes, herm_comb_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
l = legend('\Phi(k)_{Herm}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
l.Position = [0.375, 0.4, 0.1, 0.1];
legend boxoff;
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;

%ttl = title('C', 'FontSize',15, 'FontWeight', 'bold'); 
% ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
ax.FontSize = 12;
ax.LineWidth = 2;
ttl = title('C', 'FontSize',14, 'FontName','Arial','FontWeight','normal'); 
ttl.Units = 'Normalize'; ttl.Position(1) = -0.2; ttl.Position(2) = 1.05;

set(gcf, 'Color', 'W');
set(gca,'box','off');

 
%male weighted
subplot(2,2,4);
plot(male_comb_wei_rcc,'r','LineWidth',2);
hold on;
plot(male_comb_wei_mean, '-g','LineWidth',2);
plot(male_comb_norm_wei_rcc,'b','LineWidth',1.5);
plot(male_comb_wei_rich_nodes, male_comb_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
l = legend('\Phi(k)_{Male}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
l.Position = [0.695, 0.35, 0.1, 0.1];
legend boxoff;
ax= gca; 
ax.FontSize = 12;
ax.LineWidth = 2;
ttl = title('D', 'FontSize',14, 'FontName','Arial','FontWeight','normal'); 
ttl.Units = 'Normalize'; ttl.Position(1) = -0.2; ttl.Position(2) = 1.05;

set(gcf, 'Color', 'W');
set(gca,'box','off');

% Set the size of the figure
width = 900;  % Width in pixels
height = 600; % Height in pixels
set(fig, 'Position', [100, 100, width, height]);

lgd.Position = [0.2, 0.2, 0.2, 0.2];
%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "rich_club_comb_cook.pdf");
%export_fig(fig_save_path);
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');



%% chem plots

fig = figure();
%herm unweighted
subplot(2,2,1);
grey  = [220,220,220]./255;
x = [48 68 68 48];
cook_herm_chem_rich_neurons = cook_herm_chem_lf.cell_name(x(1)<=cook_herm_chem_lf.degree & cook_herm_chem_lf.degree<=x(2));

y = [0.02 0.02 2.5 2.5];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(herm_chem_bin_rcc,'r','LineWidth',2); 
plot(herm_chem_bin_mean, '-g','LineWidth',2);
plot(herm_chem_norm_bin_rcc,'-b','LineWidth',1.5);
plot(herm_chem_rich_nodes,herm_chem_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K');
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1, 'LineWidth', 1.5);
legend('\Phi(k))_{Herm}', '\Phi(k))_{random}', '\Phi(k))_{normalized}', 'Location', 'northeast');
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('A', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
%title('Unweighted rich club curve(herm)'); 


%male unweighted
subplot(2,2,2);
grey  = [220,220,220]./255;
x = [40 62 62 40];
cook_male_chem_rich_neurons = cook_male_chem_lf.cell_name(x(1)<=cook_male_chem_lf.degree & cook_male_chem_lf.degree<=x(2));

y = [0.02 0.02 3 3];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(male_chem_bin_rcc,'-r','LineWidth',2); 
plot(male_chem_bin_mean, '-g','LineWidth',2);
plot(male_chem_norm_bin_rcc,'-b','LineWidth',1.5);
plot(male_chem_rich_nodes,male_chem_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K'); 
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1,'LineWidth', 1.5);
legend('\Phi(k)_{Male}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('B', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');


%herm weighted 
subplot(2,2,3);
plot(herm_chem_wei_rcc,'r','LineWidth',2);
hold on;
plot(herm_chem_wei_mean, '-g','LineWidth',2);
plot(herm_chem_norm_wei_rcc,'b','LineWidth',1.5);
plot(herm_chem_wei_rich_nodes, herm_chem_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
legend('\Phi(k)_{Herm}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
%ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('C', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
set(gca,'box','off');

 
%male weighted
subplot(2,2,4);
plot(male_chem_wei_rcc,'r','LineWidth',2);
hold on;
plot(male_chem_wei_mean, '-g','LineWidth',2);
plot(male_chem_norm_wei_rcc,'b','LineWidth',1.5);
plot(male_chem_wei_rich_nodes, male_chem_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
legend('\Phi(k)_{Male}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('D', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
set(gca,'box','off');
%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "rich_club_chem_cook.pdf");
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');

%% elec plots
fig = figure();
%herm unweighted
subplot(2,2,1);
grey  = [220,220,220]./255;
x = [18 33 33 18];
y = [0.02 0.02 2 2];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(herm_elec_bin_rcc,'r','LineWidth',2); 
plot(herm_elec_bin_mean, '-g','LineWidth',2);
plot(herm_elec_norm_bin_rcc,'-b','LineWidth',1.5);
plot(herm_elec_rich_nodes,herm_elec_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K');
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1, 'LineWidth', 1.5);
legend('\Phi(k)_{Herm}', '\Phi(k)_{random}', '\Phi(k)_{normalized}', 'Location', 'northeast');
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('A', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
%title('Unweighted rich club curve(herm)'); 


%male unweighted
subplot(2,2,2);
grey  = [220,220,220]./255;
x = [12 42 42 12];
y = [0.02 0.02 2.5 2.5];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(male_elec_bin_rcc,'-r','LineWidth',2); 
plot(male_elec_bin_mean, '-g','LineWidth',2);
plot(male_elec_norm_bin_rcc,'-b','LineWidth',1.5);
plot(male_elec_rich_nodes,male_elec_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K'); 
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1,'LineWidth', 1.5);
legend('\Phi(k)_{Male}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('B', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');


%herm weighted 
subplot(2,2,3);
plot(herm_elec_wei_rcc,'r','LineWidth',2);
hold on;
plot(herm_elec_wei_mean, '-g','LineWidth',2);
plot(herm_elec_norm_wei_rcc,'b','LineWidth',1.5);
plot(herm_elec_wei_rich_nodes, herm_elec_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
legend('\Phi(k)_{Herm}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('C', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
set(gca,'box','off');

 
%male weighted
subplot(2,2,4);
plot(male_elec_wei_rcc,'r','LineWidth',2);
hold on;
plot(male_elec_wei_mean, '-g','LineWidth',2);
plot(male_elec_norm_wei_rcc,'b','LineWidth',1.5);
plot(male_elec_wei_rich_nodes, male_elec_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
legend('\Phi(k)_{Male}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('D', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
set(gca,'box','off');
%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "rich_club_elec_cook.pdf");
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');
%% Wit 7 and wit8 plots
fig = figure();
%wit7 unweighted rc
subplot(2,2,1);
grey  = [220,220,220]./255;
x = [26 31 31 26];
wit7_rich_neurons = wit7_lf.cell_name(x(1)<=wit7_lf.degree & wit7_lf.degree<=x(2));
wit7_rich_cell_class = wit7_lf.cell_class(x(1)<=wit7_lf.degree & wit7_lf.degree<=x(2));

y = [0.02 0.02 2 2];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(wit7_bin_rcc,'r','LineWidth',2); 
plot(wit7_bin_mean, '-g','LineWidth',2);
plot(wit7_norm_bin_rcc,'-b','LineWidth',1.5);
plot(wit7_rich_nodes,wit7_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K');
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1, 'LineWidth', 1.5);
legend('\Phi(k)_{Wit7}', '\Phi(k)_{random}', '\Phi(k)_{normalized}', 'Location', 'northeast');
ax= gca;
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('A', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
%title('Unweighted rich club curve(herm)'); 


%wit8 unweighted
subplot(2,2,2);
grey  = [220,220,220]./255;
x = [27 38 38 27];
wit8_rich_neurons = wit8_lf.cell_name(x(1)<=wit8_lf.degree & wit8_lf.degree<=x(2));
wit8_rich_cell_class = wit8_lf.cell_class(x(1)<=wit8_lf.degree & wit8_lf.degree<=x(2));

y = [0.02 0.02 2.5 2.5];
patch(x,y,grey,'HandleVisibility','off','FaceAlpha',1, 'EdgeColor', [1,1,1]);
hold on;
plot(wit8_bin_rcc,'-r','LineWidth',2); 
plot(wit8_bin_mean, '-g','LineWidth',2);
plot(wit8_norm_bin_rcc,'-b','LineWidth',1.5);
plot(wit8_rich_nodes,wit8_rich_nodes_norm_rcc, 'ob','HandleVisibility','off'); 

xlabel('K'); 
ylabel('Unweighted rich club coefficient, \Phi');
%set(gca,'FontWeight','bold');
yline(1,'LineWidth', 1.5);
legend('\Phi(k)_{Wit8}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('B', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');


%wit7 weighted
subplot(2,2,3);
plot(wit7_wei_rcc,'r','LineWidth',2);
hold on;
plot(wit7_wei_mean, '-g','LineWidth',2);
plot(wit7_norm_wei_rcc,'b','LineWidth',1.5);
plot(wit7_wei_rich_nodes, wit7_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
legend('\Phi(k)_{Wit7}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('C', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
set(gca,'box','off');

 
%wit8 weighted
subplot(2,2,4);
plot(wit8_wei_rcc,'r','LineWidth',2);
hold on;
plot(wit8_wei_mean, '-g','LineWidth',2);
plot(wit8_norm_wei_rcc,'b','LineWidth',1.5);
plot(wit8_wei_rich_nodes, wit8_wei_rich_nodes_norm_rcc, 'ob',...
     'HandleVisibility','off');
yline(1,'LineWidth', 1.5);

xlabel('K');
ylabel('Weighted rich club coefficient, \Phi_{weighted}');
%set(gca,'FontWeight','bold');
legend('\Phi(k)_{Wit8}', '\Phi(k)_{random}', '\Phi(k)_{normalized}');
ax= gca; 
ax.FontName = 'Arial';
%ax.FontWeight = 'Bold';
ax.FontSize = 14;
ax.LineWidth = 2;
%ttl = title('D', 'FontSize',15, 'FontWeight', 'bold'); ttl.Units = 'Normalize'; ttl.Position(1) = -0.1;
set(gcf, 'Color', 'W');
set(gca,'box','off');

%%

fig_save_path = fullfile(fileparts(pwd), "figure_components", "rich_club_wit.pdf");
exportgraphics(fig,fig_save_path,'Resolution',300, 'ContentType','vector');




 