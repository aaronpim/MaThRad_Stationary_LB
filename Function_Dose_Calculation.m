function [Du] = Function_Dose_Calculation(u)
%##########################################################################
% Inputs
%   u       An [N,N,N,26,En] double tensor that represents the particle
%           probability distribution in terms of position, trajectory 
%           and energy.
% Outputs   
%   Du      An [N,N,N] tensor which denotes the dose delivered at a given 
%           point.
%##########################################################################

% Define the energy levels
    [Energy_vec] = Function_Energy_vec(size(u,5));

% Integrate over the energy variables and trajectories.
    Du = sum(trapz(Energy_vec,u,5),4);

% It is assumed that the stopping powers are identically 1.
end