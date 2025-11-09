function feat_outliers = dimorphic_neurons(feat_table, feat1, feat2)
% input: feature table and the features 
% output: neurons with 3 MAD (their name, feature values and difference

feat_diff = feat_table{:, feat1} - feat_table{:,feat2};
feat_outlier = isoutlier(feat_diff, 'median');
feat_outlier = feat_table(feat_outlier,["cell_name",feat1,feat2]);
feat_outlier.diff = feat_outlier{:,feat1}-feat_outlier{:, feat2};
feat_outliers = sortrows(feat_outlier, "diff", "descend");
end