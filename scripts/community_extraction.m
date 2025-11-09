function [NC, bestQ] = community_extraction(AM_table, rand_n)

% Input: AM: Adjacency matrix table where the first columns is neurons name
%        rand_n: number of random initializations to try
% Output: NC: a table of neurons and their community affiliation
%         bestQ: the largest modularity index
% 
neurons = AM_table{:,1};
AM = AM_table{:,2:end};

if nargin<2
    rand_n = 1000;
end

bestQ = -1;
Q_list = -1*ones(1,rand_n);

for i=2:rand_n
    rng(i)
    Mh  = 1:length(neurons);

    if i==2
        Mh  = 1:length(neurons);    % initial community affiliations
        Q_list(i)=0;            % initialize modularity value for current i
    end

    [Mh, Qi] = community_louvain(AM, [], Mh); %make just one run for community detection

    while Qi-Q_list(i-1)>1e-8          % while modularity increases
       Q_list(i-1) = Qi;
       [Mh, Qi] = community_louvain(AM, [], Mh);
    end

    Q_list(i) = Qi; %update the current iteration modularity index



    if Qi > bestQ 
        bestQ = Qi;
        best_Mh = Mh;
    end
end

NC = table(neurons, best_Mh, 'VariableNames', {'Neuron', 'Community'});

end