% Generates a 9x9 correlation matrix from the original 83x83 Empirical
% correlation matrix. By grouping the 83 regions according to the larger
% networks they are a part of (ie. Visual Network, DMN, etc.) and
% averaging.

clear; clc;

load regions_83_to_9.mat % Lists indices for 9 networks
load Corr_FMRI.mat % Empirical correlation matrix

row = []; % Contains all the indices of a particular region
Corr_FMRI_9FN = zeros(size(region,1));

for i = 1:size(region,1)
    
    row = region(i,:);
    Lrow = logical(row);
    row(~Lrow) = [];
    
    for j = i:size(region,1)
        
        col = region(j,:);
        Lcol = logical(col);
        col(~Lcol) = [];
        
        Corr_FMRI_9FN(i,j) = mean(mean(Corr_FMRI(row,col)));
        
    end

end

Corr_FMRI_9FN = triu(Corr_FMRI_9FN) + triu(Corr_FMRI_9FN,1)';
imagesc(Corr_FMRI_9FN)

clear col Lcol row Lrow i j 