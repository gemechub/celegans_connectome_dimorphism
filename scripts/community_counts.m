function [unique_communities, unique_features, feature_counts] = community_counts(NC, anatomical_info, feature)

[commonNeurons, idx_NC, idx_anatomical] = intersect(NC.Neuron, anatomical_info.cell_name, 'stable');

community = NC.Community(idx_NC);
features = anatomical_info{idx_anatomical,feature};

unique_communities = unique(community);
unique_features = unique(features);

if strcmp(feature, 'cell_type') && length(unique_features)==3
    unique_features = {'sensory', 'interneuron', 'motor'}; %just to have this order
elseif strcmp(feature, 'cell_type') && length(unique_features)==4 %this is from Witv dataset
    unique_features = {'Sensory neuron', 'Interneuron', 'Modulatory neuron','Motor neuron'};
elseif strcmp(feature, 'general_location')
    unique_features = {'tail', 'mid-body', 'head'};

end

feature_counts = zeros(numel(unique_communities), numel(unique_features));
% Count occurrences of each general_location within each community
for i = 1:numel(unique_communities)
    for j = 1:numel(unique_features)
        % Find neurons in current community and location
        feature_counts(i, j) = sum(community == unique_communities(i) & strcmp(features, unique_features{j}));
    end
end

end