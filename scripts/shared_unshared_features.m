function [shared_neur_feats, unshared_neur_feats] = shared_unshared_features(features_table, feature, sex)
% Input: feature_table that can be accessed with .feature 
%        feature: feature (eg: degree) for which we want to extract feature values
% Output: shared_neur_feats: features of neurons shared by both sexes
%         unshared_neur_feats: features of neurons that sex specific
if strcmp(sex, 'male')
   sex_specific_neurons = strcmp(features_table.Gender, 'MS');
elseif strcmp(sex,'herm')
    sex_specific_neurons = strcmp(features_table.Gender, 'HS');
end
sex_shared_neurons = strcmp(features_table.Gender, 'GS');
%feature_list = features_table(:,feature);
%shared_neur_feats = feature_list(sex_shared_neurons);
%unshared_neur_feats = feature_list(sex_specific_neurons);
shared_neur_feats = features_table(sex_shared_neurons, feature);
unshared_neur_feats = features_table(sex_specific_neurons, feature);
end