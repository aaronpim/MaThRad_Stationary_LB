function [data] = Dual_Data_Calculation(u,k,S,alpha,u_t)
%##########################################################################
% Inputs
%   u       An [N,N,N,26,En] double tensor that represents the particle
%           probability distribution in terms of position, trajectory 
%           and energy.
%   k       An [N,N,N,En] double tensor that represents the stopping powers
%           of a particle in a given medium and for a given energy.
%   S       An [N,N,N,M] logical tensor that represents the spatial
%           position of the M subdomains
%   alpha   An [M,1] double vector that represents the weights in each
%           individual region
%   u_t     An [N,N,N,M] double tensor that represents the target dose in
%           each of the subregions.
% Outputs   
%   data    An [N,N,N,26,En] tensor which denotes the dose delivered at a  
%           given point.
%##########################################################################
if nargin == 1
    k = ones(size(u,[1,1,1,5]));
end
% Calculate the dose
    [Du] = Function_Dose_Calculation(u,k);

% Calculate the data for the dual equation.
data = 0;
for i = 1:length(alpha)
    data = data + alpha(i)*S(:,:,:,i)*(Du-u_t(:,:,:,i));
end
    data = k.*repmat(data,[1,1,1,En]);
    data = repmat(data,[1,1,1,1,size(u,4)]);
    data = permute(data,[1,2,3,5,4]);
end