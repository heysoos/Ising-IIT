% Generate a 6x6 matrix from an 83x83 size by grouping regions of the 83x83
% matrix as defined by indice_map.mat which is a matrix of n*m (6x? at the
% moment). The nodes in each region are then averaged over their connection
% with the other regions.

% This code was used to compress a large (83x83) empirical correlation of
% the brain (averaged over (17?) subjects) into a smaller (6x6) matrix.
% This was done to generate a functional connectivity matrix and use it as
% an parameter for an Ising model simulation.



clear all
clc

% the matrix representing the regions in the brain and its nodes
load indice_map.mat
% the emprical correlation matrix of the brain (83x83)
load MJ_all.mat

J_brain = zeros(6);

for i = 1:size(indice_map,1)
    ind1 = indice_map(i,:);
    
    for j = (i+1):size(indice_map,1);
        ind2 = indice_map(j,:);
        J_brain(i,j) = mean(mean(MJ_all_mn(ind1,ind2)));
        % averages the connectivity between two networks
    end
end

J_brain = triu(J_brain) + triu(J_brain,1)';

% Each row represents the following regions: Auditory, DMN, ECN (left and
% right), Salience, Sensorimotor, Visual (L,M,O) for a total of 6 regions.

% NOTE: These regions have overlapping indices. For example the DMN, ECN(l
% and r) and Visual (L,M, and O) regions both have indices 20 and 61 in
% their regions. Regions 20 and 61 are the PRECUNEUS left and right
% respectively.