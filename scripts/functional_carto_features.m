function com_features = functional_carto_features(AM_table, NC, cell_features)
%input : AM_table : adjacency matrix table where the first column is name of the neurons
%       NC : a table containing all cells and their community affiliaton
%       cell_features: a table containing all cells and their basic metafeaures such as cell type
%
%output : a table containing the local features of each neuron in the AM_table

    %extract name of neurons and convert the AM_table to AM
    neurons = AM_table{:,1};
    AM = AM_table{:,2:end};
    com_vector = NC.Community;

    % Sort the com_vector to be similar to AM neurons
    % Reorder com_vector based on neurons list
    [~, idx] = ismember(neurons, NC.Neuron); % Find indices of neurons in com_vector
    sorted_com_vector = com_vector(idx, :);              % Reorder rows of com_vector

    % Get the neurons features
    [~,idx] = ismember(neurons, cell_features.cell_name);
    cell_type = cell_features.cell_type(idx,:);
    cell_sex = cell_features.Gender(idx,:);

    %% calculate the within-module degree z-score and participation coefficient
    module_zscore = module_degree_zscore(AM, sorted_com_vector, 3);
    participation_coefficient = participation_coef(AM, sorted_com_vector, 1);

    %% For each cell, assign region in the functional cartography 
    region = zeros(length(module_zscore),1);
    for i=1:length(module_zscore)
        modz = module_zscore(i);
        partcoef = participation_coefficient(i);
        if partcoef<0.05 && modz<2.5
            r = 1;
        elseif 0.05<=partcoef && partcoef<0.625 && modz<2.5
            r = 2;
        elseif 0.625<= partcoef && partcoef<0.8 && modz<2.5
            r = 3;
        elseif partcoef>=0.8 && modz<2.5
            r = 4;
        elseif partcoef<0.3 && modz>=2.5
            r = 5;
        elseif 0.3<=partcoef && partcoef<0.75 && modz>=2.5
            r = 6;
        elseif partcoef>=0.75 && modz>=2.5
            r = 7;
        end

        region(i)=r;
    end


    %% make the final table
    com_features = table(neurons, cell_type, cell_sex, module_zscore, participation_coefficient, region);
  
end
