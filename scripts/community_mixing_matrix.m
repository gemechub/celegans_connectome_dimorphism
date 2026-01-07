function [mixTbl, mixCountTbl,mixRowNormTbl,mixGlobalNormTbl,edgeList, K ] = community_mixing_matrix(commTbl, edgeTbl)

%% Inputs
%   commTbl : community table that has a neuron (Neuron column) and its community assign in 'Community' column
%   edgeTbl : the edge list table that 'Source', 'Target' and 'Weight' columns
% Outputs
%   mixTbl : mixing matrix where each cell is the sum of weights of edges between those two 
%           communities. Direction is from row to column
%   mixCountTbl : mixing matrix where each cell is the count of edges between those two communities
%   mixRowNormTbl : mixTbl, but normalized across rows
%   mixGlobalNormTbl : mixTbl, but normalized by sum of all weights
%   edgeList : edgelist version of mixTbl
%   K : community names


% Column names
neuronCol    = "Neuron";
communityCol = "Community";

srcCol   = "Source";
tgtCol   = "Target";
wCol     = "Weight";

%% Ensure neuron identifiers are strings for consistent matching
commTbl.(neuronCol) = string(commTbl.(neuronCol));
edgeTbl.(srcCol)    = string(edgeTbl.(srcCol));
edgeTbl.(tgtCol)    = string(edgeTbl.(tgtCol));

%% Build a map: neuron name -> community label (can be numeric or string)
% Normalize communities to categorical indices 1..K
commLabels = commTbl.(communityCol);

if isnumeric(commLabels) || islogical(commLabels)
    commLabelsCat = categorical(commLabels);
else
    commLabelsCat = categorical(string(commLabels));
end

% Community names (for labeling M)
commNames = categories(commLabelsCat);
K = numel(commNames);

% Map each neuron to its community index (1..K)
neuronToCommIdx = containers.Map( ...
    cellstr(commTbl.(neuronCol)), ...
    num2cell(double(commLabelsCat)) ...
);

%% For each edge, look up source/target community indices
nE = height(edgeTbl);

srcComm = nan(nE, 1);
tgtComm = nan(nE, 1);

for e = 1:nE
    s = edgeTbl.(srcCol)(e);
    t = edgeTbl.(tgtCol)(e);

    if isKey(neuronToCommIdx, char(s))
        srcComm(e) = neuronToCommIdx(char(s));
    end
    if isKey(neuronToCommIdx, char(t))
        tgtComm(e) = neuronToCommIdx(char(t));
    end
end

% Keep only edges where both endpoints have community labels
keep = ~isnan(srcComm) & ~isnan(tgtComm);

srcComm = srcComm(keep);
tgtComm = tgtComm(keep);

% Weights: if missing, default to 1 (edge count behavior)
if ismember(wCol, edgeTbl.Properties.VariableNames)
    w = edgeTbl.(wCol)(keep);
    w = double(w);
else
    w = ones(sum(keep), 1);
end

%% Construct directed weighted mixing matrix M (K x K)
% M(a,b) = sum of weights for edges from community a -> community b
M = accumarray([srcComm, tgtComm], w, [K, K], @sum, 0);

% Also build edge-count mixing matrix (unweighted)
Mcount = accumarray([srcComm, tgtComm], 1, [K, K], @sum, 0);

%% Optional normalizations
% Row-normalized: for each source community a, proportions of outgoing weight to target communities
rowSums = sum(M, 2);
M_rowNorm = M ./ max(rowSums, eps);

% Global-normalized: fraction of total weight
M_globalNorm = M ./ max(sum(M, "all"), eps);

%% Wrap results in a table for readability
mixTbl = array2table(M, "VariableNames", commNames, ...
                        "RowNames", commNames);

mixCountTbl = array2table(Mcount, "VariableNames", commNames, ...
                               "RowNames", commNames);

mixRowNormTbl = array2table(M_rowNorm, "VariableNames", commNames, ...
                                  "RowNames", commNames);

mixGlobalNormTbl = array2table(M_globalNorm, "VariableNames", commNames, ...
                                     "RowNames", commNames);

%% Get edgelist version of the tables
%% Community names: Comm1, Comm2, ..., CommK
commNames = "Comm" + string(1:K);

%% Build edge list
[srcIdx, tgtIdx, w] = find(M);   % directed edges with nonzero weight

edgeList = table( ...
    commNames(srcIdx)', ...
    commNames(tgtIdx)', ...
    w, ...
    'VariableNames', {'Source','Target','Weight'} ...
);


end