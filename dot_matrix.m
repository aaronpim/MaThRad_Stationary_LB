function [vals,assign_matrix] = dot_matrix(s_vec)
%##########################################################################
% Inputs
%   s_vec           An [M,3] double matrix which denotes the possible beam
%                   trajectories.
% Outputs
%   vals            A vector which contains all of the unique dot products 
%                   of the normalised s_vec.
%   assign_matrix   An [M,M] matrix of integers which indexes which dot
%                   product value in vals each pair of s_vec rows 
%                   corresponds to.
%##########################################################################

% Normalise the trajectory vectors
    S = normr(s_vec);

% Compute all possible dot products
    dot_matx = S*S';

% Select unique values
% Accounting for machine error
    vals = uniquetol(S*S',0.001);

% Define assign_matrix
    assign_matrix = zeros(size(dot_matx));
    for i = 1:length(vals)
        assign_matrix = assign_matrix + i*(abs(dot_matx - vals(i)) <0.001);
    end
end