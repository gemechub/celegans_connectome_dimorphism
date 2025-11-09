function [similarity_score, same_diff_matrix1, same_diff_matrix2] = community_similarity(NC1, NC2)
% Inputs
% Outputs

% Step 1: Find overlapping neurons
common_neurons = intersect(NC1.Neuron, NC2.Neuron);

% Step 2: Sort the tables based on common neuron order
[~, idx1] = ismember(common_neurons, NC1.Neuron);
[~, idx2] = ismember(common_neurons, NC2.Neuron);


communities1 = NC1.Community(idx1);
communities2 = NC2.Community(idx2);

% Number of overlapping neurons
num_neurons = length(common_neurons);
same_diff_matrix1 = zeros(num_neurons);
same_diff_matrix2 = zeros(num_neurons);

for i = 1:num_neurons
    for j = 1:num_neurons
        same_diff_matrix1(i, j) = (communities1(i) == ...
            communities1(j));

        same_diff_matrix2(i, j) = (communities2(i) == ...
            communities2(j));
    end
end


non_zero_pairs = (same_diff_matrix1 + same_diff_matrix2) > 0;

%cook_herm_male
matching_entries = sum(sum((same_diff_matrix1 == ...
    same_diff_matrix2) & non_zero_pairs));
total_non_zero_entries = sum(non_zero_pairs(:));  % Total pairs with at least one non-zero

if total_non_zero_entries == 0
    similarity_score = 0;
else
    similarity_score = matching_entries / total_non_zero_entries;
end

end