function [fitted_ccdf, ksh, ksp,fitted_features, fitted_params] = fit_distribution(feature, dist_type)
% fit_distribution fits a given feature vector to different distributions 
% and returns the Complementary Cumulative Distribution Function (CCDF).
%
% Inputs:
%   feature    - A vector of values (e.g., degrees, strengths, or other features).
%   dist_type  - The type of distribution to fit ('logn', 'wbl', 'exp', 'power').
%
% Outputs:
%   unique_feature - Unique feature values used for fitting.
%   fitted_ccdf    - Fitted Complementary CDF values.
%   ksh            - Kolmogorov-Smirnov (KS) test result (1 = reject null, 0 = fail to reject).
%   ksp            - KS test p-value.
%   fitted_features - the min to max vector of the unique features
%   fitted_params - parameters of the fitted distribution


% Remove non-positive values (valid filtering for all fits)**
valid_feature = feature(feature > 0);
if isempty(valid_feature)
    error('No valid data points remain after filtering non-positive values.');
end

% Extract unique sorted values
unique_feature = unique(valid_feature); 

fitted_features = 0:max(unique_feature);

% Fit the specified distribution
alpha = NaN;  % Default NaN (only updated for power-law)

if strcmp(dist_type, 'logn')
    % Lognormal distribution fitting
    fit_params = lognfit(valid_feature);
    mu = fit_params(1);
    sigma = fit_params(2);
    %cdf = logncdf(min(unique_feature):max(unique_feature),mu,sigma);
    cdf = logncdf(fitted_features,mu,sigma);
    feature_pd = makedist('Lognormal', 'mu', mu, 'sigma', sigma);
    fitted_params.mu = mu;
    fitted_params.sigma = sigma;

elseif strcmp(dist_type, 'wbl')
    % Weibull distribution fitting
    fit_params = wblfit(valid_feature);
    a = fit_params(1);
    b = fit_params(2);
    cdf = wblcdf(fitted_features,a,b);
    feature_pd = makedist('Weibull', 'a', a, 'b', b);
    fitted_params.a = a;
    fitted_params.b = b;
    

elseif strcmp(dist_type, 'exp')
    % Exponential distribution fitting
    mu = expfit(valid_feature);  % Use mean as the parameter
    cdf = expcdf(fitted_features,mu);
    feature_pd = makedist('Exponential', 'mu', mu);
    fitted_params.mu = mu;

elseif strcmp(dist_type, 'power')
    % Estimate xmin (minimum value for power-law behavior)
    xmin = min(valid_feature); % Lowest value in valid_feature

    % Compute the MLE of alpha
    n = length(valid_feature);
    alpha = 1 + n / sum(log(valid_feature / xmin)); % MLE estimate

    % Compute the fitted CDF
    cdf = (min(unique_feature):max(unique_feature) / xmin).^(- (alpha - 1));  % Power-law CCDF formula

    % Define a custom CDF function for KS test
    feature_pd = @(x) (x / xmin).^(- (alpha - 1));

    fitted_params.alpha = alpha;

else
    error('Unsupported distribution type. Use "logn", "wbl", "power", or "exp".');
end


% Perform KS test
if strcmp(dist_type, 'power')
    % Compute the fitted CCDF
    fitted_ccdf = cdf;
    fitted_features = min(unique_feature):max(unique_feature);
    test_ccdf = (valid_feature / xmin).^(- (alpha - 1));
    ks_test_values = [valid_feature, 1-test_ccdf];
    [ksh, ksp,kss] = kstest(valid_feature, 'CDF',ks_test_values);
else % test using probability distribution object
    % Compute the fitted CCDF
    fitted_ccdf = 1 - cdf;
    [ksh, ksp, kss] = kstest(valid_feature, 'CDF', feature_pd);
end



end
