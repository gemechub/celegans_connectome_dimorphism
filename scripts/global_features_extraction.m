function [global_feat_table_binary, global_feat_table_weighted] = global_features_extraction(AM, table)
% input: AM : adjacency matrix array or table where the first column is name of the neurons
%       network_type: 'binary' or 'weighted'. 
%       table: some of the AMs are in table format while some are in array format, 
%               if it is in table format set table==1, if array 0
% output: global_features_table

%convert the AM_table to AM
    if table==1
        AM = AM{:,2:end};
    end

    % Binary global features
    no_of_features = 8;
    %binarize the adjacency matrix (1 if there is a connection 0 if not)
    BM = weight_conversion(AM, 'binarize');
    %get the distance matrix
    DM = distance_bin(AM);
    
    glob_features = zeros(1, no_of_features);
    
    for j=1:no_of_features
        switch j
            case 1 %assortavity
                glob_features(j) = assortativity_bin(BM, 0);
            case 2 %cpl
                glob_features(j) = charpath(DM);
            case 3 %radius
                [~,~,~, radius, ~] = charpath(DM);
                glob_features(j) = radius;
            case 4 %diameter
                [~,~,~,~,diameter] = charpath(DM);
                glob_features(j) = diameter;
            case 5 %transitivity
                glob_features(j) = transitivity_bd(BM);
            case 6 %density
                glob_features(j) = density_dir(BM);
            case 7 %global efficiency
                glob_features(j) = efficiency_bin(BM);
            case 8 %mean clustering coefficient
                cc = clustering_coef_bd(BM);
                glob_features(j) = mean(cc(~isinf(cc)));
        end
        
    end
    
    %convert it to table with features as column header
    features = {'assortavity', 'cpl', 'radius', 'diameter',...
        'transitivity', 'density','global_efficiency', 'mean_cc'};
    global_feat_table_binary = array2table(glob_features, 'VariableNames', features);
   

    %weighted local features extraction
    no_of_features = 8;
    %get length matrix, which is ued for some measures
    LM = weight_conversion(AM, 'lengths');
    %get the distance matrix
    DM = distance_wei(LM);
    %normalize for efficiency measure
    NM = weight_conversion(AM, 'normalize');
    
    glob_features = zeros(1, no_of_features);

    for j=1:no_of_features
        switch j
            case 1 %assortavity
                glob_features(j) = assortativity_wei(AM, 0);
            case 2 %cpl
                glob_features(j) = charpath(DM);
            case 3 %radius
                [~,~,~, radius, ~] = charpath(DM);
                glob_features(j) = radius;
            case 4 %diameter
                [~,~,~,~,diameter] = charpath(DM);
                glob_features(j) = diameter;
            case 5 %transitivity
                glob_features(j) = transitivity_wd(NM);
            case 6 %density
                glob_features(j) = density_dir(AM);
            case 7 %global efficiency
                glob_features(j) = efficiency_wei(AM);
            case 8 %mean clustering coefficient
                glob_features(j) = mean(clustering_coef_wd(AM));
        end
    end
    %convert it to table with features as column header
    features = {'assortavity', 'cpl', 'radius', 'diameter',...
    'transitivity', 'density','global_efficiency', 'mean_cc'};
    global_feat_table_weighted = array2table(glob_features, 'VariableNames', features);
        

end

