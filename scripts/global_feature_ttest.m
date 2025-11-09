function [h, p] = global_feature_ttest(null_models_gf, original_gf, feature_names)
%% Input: null_models_gf - table of null models global features
%         original_gf - the table containing the original connectome's global feature values
%                       to be compared against nulls
%         feature_names - the global features' names
% Output: h - the null hypothesis test result: 1 is rejected, 0 accepted
%         p - p-value of the test
for i=1:length(feature_names)
    feature_name = feature_names{i};
    feature_null_values = null_models_gf{:,feature_name};
    feature_null_values = feature_null_values(~isnan(feature_null_values));
    feature_null_values = feature_null_values(~isinf(feature_null_values));
    [h,p] = ttest(feature_null_values, original_gf{:,feature_name});
    fprintf('%s: %d\n', feature_name, p);
 
end
