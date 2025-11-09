function local_feat_table = local_features_extraction(AM_table, cell_type_table, column_name)
%input : AM_table : adjacency matrix table where the first column is name of the neurons
%        cell_type_table : a table containing all cells and their basic meta features
%        column_name : either 'cell_name' or 'cell_class' because Cook and Witvliet neurons names are a bit different
%
%output : a table containing the local features of each neuron in the AM_table

    %extract name of neurons and convert the AM_table to AM
    neurons = AM_table{:,1};
    AM = AM_table{:,2:end};

    %get the meta neurons features
    if nargin < 3
        column_name = 'cell_name';
    end

    if strcmp(column_name, 'cell_name')
        [~,loc_of_neurons] = ismember(neurons,cell_type_table.cell_name);
        meta_features = cell_type_table(loc_of_neurons,:);
    elseif strcmp(column_name,'cell_class')
        neurons_class = cellfun(@(x) x(1:3), neurons, 'UniformOutput', false);
        [~,loc_of_neurons] = ismember(neurons_class,cell_type_table.cell_class);
        meta_features = cell_type_table(loc_of_neurons,:);
        meta_features.cell_name = neurons;
    end

    

    %get length matrix, which is ued for some measures
    no_of_nodes = length(AM);
    herm_LM = weight_conversion(AM, 'lengths');
    %get the distance matrix
    herm_DM = distance_wei(herm_LM);
    %normalize for efficiency measure
    herm_NM = weight_conversion(AM, 'normalize');
    herm_BM = weight_conversion(AM, 'binarize');
    
    no_of_features = 15;

    loc_features = zeros(no_of_nodes, no_of_features);
            
    %extract the the features and save it in a table
    for j=1:no_of_features
        switch j
            case 1 %indegree
                [indegree,~,~] = degrees_dir(AM);
                loc_features(:,j) = indegree';
            case 2 %outdegree
                [~,outdegree,~] = degrees_dir(AM);
                loc_features(:,j) = outdegree';
            case 3 %degree
                [~,~,degree] = degrees_dir(AM);
                loc_features(:,j) = degree';
            case 4 %in strength
                [instrength,~,~] = strengths_dir(AM);
                loc_features(:,j) = instrength';
            case 5 %outstrength
                [~,outstrength,~] = strengths_dir(AM);
                loc_features(:,j) = outstrength';
            case 6 %total strength
                [~,~,strength] = strengths_dir(AM);
                loc_features(:,j) = strength';
            case 7 %betweenness centrality
                bc_wei = betweenness_wei(herm_LM);
                loc_features(:,j) = bc_wei;
            case 8 %weighted betweenness centrality normalized
                bc_wei_norm = bc_wei/((no_of_nodes-1)*(no_of_nodes-2));
                loc_features(:,j) = bc_wei_norm;
            case 9 % betweenness centrality binary
                bc_bin = betweenness_bin(herm_BM);
                loc_features(:,j) = bc_bin;
            case 10 %betweenness centrality normalized 
                bc_bin_norm = bc_bin/((no_of_nodes-1)*(no_of_nodes-2));
                loc_features(:,j) = bc_bin_norm;
            case 11 %eccentrcity
                [~,~,curr_subj_ecc,~,~] = charpath(herm_DM); 
                loc_features(:,j) = curr_subj_ecc;
            case 12 %clustering coefficient
                loc_features(:,j) = clustering_coef_wd(herm_NM);
            case 13 %clustering coefficient binary
                loc_features(:,j) = clustering_coef_bd(herm_BM);
            case 14 %local effeciecny
                loc_features(:,j) = efficiency_wei(herm_NM, 2);
            case 15 %local efficiency binay
                loc_features(:,j) = efficiency_bin(herm_BM,1);
        end
    end
        %convert it to table with features as column header
    features = {'indegree', 'outdegree', 'degree', 'instrength',... 
        'outstrength', 'strength', 'BC_wei', 'BC_wei_norm','BC_bin',...
        'BC_bin_norm', 'eccentricity','CC_wei', 'CC_bin', 'LE_wei', 'LE_bin'};
    local_feat_table = array2table(loc_features, 'VariableNames', features);
    local_feat_table = [meta_features, local_feat_table];
end
