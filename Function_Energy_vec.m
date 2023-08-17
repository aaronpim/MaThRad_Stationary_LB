function [Energy_vec] = Function_Energy_vec(En)
%##########################################################################
% Inputs
%   En          An integer value which denotes the number of energy levels.
% Outputs   
%   Energy_vec  An [En,1] double vector which indexes the energy levels.
%##########################################################################

% Define Energy_vec
%Energy_vec = 10.^(linspace(-4,log10(350/469),En));
Energy_vec = linspace(1.0e-4,350/469,En);

% This vector is used in multiple points and thus to easy editing, I
% define this as an overhead variable. This assumes that the maximum energy
% a proton can have is 350 MeV.
end