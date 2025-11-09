function region_count = functional_carto_regions_count(modz, partcoef)
% counts the number of neurons in each region of the functional cartography
    r1 = nnz(partcoef<0.05 & modz<2.5);
    r2 = nnz(0.05<=partcoef & partcoef<0.625 & modz<2.5);
    r3 = nnz(0.625<= partcoef & partcoef<0.8 & modz<2.5);
    r4 = nnz(partcoef>=0.8 & modz<2.5);
    r5 = nnz(partcoef<0.3 & modz>=2.5);
    r6 = nnz(0.3<=partcoef & partcoef<0.75 & modz>=2.5);
    r7 = nnz(partcoef>=0.75 & modz>=2.5);
    region_count = table(r1,r2,r3,r4,r5,r6,r7);
end
