function [sigma_a] = Function_sigma_a_dist(Inclusion)
%##########################################################################
% Inputs
%   Inclusion   An [N,N,N,26,En] logical tensor which denotes the position
%               of the dense inclusion.
% Outputs
%   sigma_a     An [N,N,N,26,En] double tensor which denotes the absorpsion
%               coefficient at each point, each trajectory and each
%               energy.
%##########################################################################

% Define number of energy levels
En = size(Inclusion,5);

% Define the energy vec
[Energy_vec] = Function_Energy_vec(En);

% Define sigma_a in the inclusion
Sa_inc = @(i) 0.0*(1-Energy_vec(i));

% Define sigma_a outside of the inclusion
Sa_exc = @(i) 0.0*(1-Energy_vec(i));

% Pre-allocate sigma_a
sigma_a = zeros(size(Inclusion));

% Define sigma_a(.,E)
for E = 1:size(Inclusion, 5)
    sigma_a(:,:,:,:,E) = Sa_inc(E)*Inclusion(:,:,:,:,E)+...
                            Sa_exc(E)*(~Inclusion(:,:,:,:,E));
end
end