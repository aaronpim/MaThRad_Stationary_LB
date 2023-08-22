function [Du] = Function_Dose_Calculation(u,k)
%##########################################################################
% Inputs
%   u       An [N,N,N,26,En] double tensor that represents the particle
%           probability distribution in terms of position, trajectory 
%           and energy.
%   k       An [N,N,N,En] double tensor that represents the stopping powers
%           of a particle in a given medium and for a given energy.
% Outputs   
%   Du      An [N,N,N] tensor which denotes the dose delivered at a given 
%           point.
%##########################################################################
if nargin == 1
    k = ones(size(u,[1,1,1,5]));
end
% Define the energy levels
    [Energy_vec] = Function_Energy_vec(size(u,5));

% Integrate over the trajectory variable
% (4*pi^2/(2+3*size(u,4))) = area of each face on the surface of a sphere.
    u_temp = (4*pi^2/(2+3*size(u,4)))*reshape(sum(u,4),size(u,[1,2,3,5]));
    
% Integrate over the energy variable
    Du = trapz(Energy_vec,k.*u_temp,4);
end