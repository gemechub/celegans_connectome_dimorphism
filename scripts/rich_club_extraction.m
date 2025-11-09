function [bin_rcc, bin_rcc_null, norm_bin_rcc, wei_rcc, wei_rcc_null, norm_wei_rcc] = ...
    rich_club_extraction(AM_table, null_AMs)

% inputs: AM_table: the original connectomes AM table
%       : null_AMs: a 3D array containing the null models
% outputs: bin_rcc: binary rich club coefficients for the original AM
%        : bin_rcc_null: binary rich club coefficients for the null models
%        : norm_bin_rcc: normalized binary rich club coefficients
%        : wei_rcc: weighted rich club coefficients
%        : wei_rcc_null: weighted rich club coefficients for the null models
%        : norm_wei_rcc: normalized weighted rich club coefficients

AM = AM_table{:, 2:end}; 

bin_rcc_null = [];
wei_rcc_null = [];

for i=1:1000
    if mod(i,100)==0
        fprintf('%d... of 1000\n', i);
    end
    bin_rcc_null(i,:) = rich_club_bd(weight_conversion(squeeze(null_AMs(i,:,:)), 'binarize'));
    wei_rcc_null(i,:)= rich_club_wd(squeeze(null_AMs(i,:,:)));
end
 
 bin_rcc = rich_club_bd(weight_conversion(AM, 'binarize'));
 wei_rcc = rich_club_wd(AM);
 
 
 bin_rcc_null_mean = mean(bin_rcc_null, 'omitnan');
 wei_rcc_null_mean = mean(wei_rcc_null, 'omitnan');
 
 norm_bin_rcc = bin_rcc./bin_rcc_null_mean;
 norm_wei_rcc = wei_rcc./wei_rcc_null_mean;