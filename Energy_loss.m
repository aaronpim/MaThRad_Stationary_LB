function [output] = Energy_loss(u, S, adjointyn)
%##########################################################################
% Inputs
%   S           An [N,N,N,M,En] double tensor which denotes the mean energy
%               loss per unit length.
%   u           An [N,N,N,26,En] double tensor that represents the particle
%               probability distribution in terms of position, trajectory 
%               and energy.
%   adjointyn   A logical value to determine if we consider the operator or
%               it's adjoint.
% Outputs
%   outputs     An [N,N,N,M,En] double tensor which represents the change
%               in energy operator applied to u.
%##########################################################################

% Default consider the forward operator.
if nargin == 2
    adjointyn = false;
end

% Precompute the output
output = zeros(size(u)); En = size(u,5);

% Define the energy levels
[Energy_vec] = Function_Energy_vec(En);

if ~adjointyn
    for E = 1:En-1
        output(:,:,:,:,E) = output(:,:,:,:,E) + (S(E+1)*u(:,:,:,:,E+1) - S(E)*u(:,:,:,:,E))./(Energy_vec(E+1)-Energy_vec(E));
    end
else
    for E = 1:En-1
        output(:,:,:,:,E) = output(:,:,:,:,E) + S(E)*(u(:,:,:,:,E+1) - u(:,:,:,:,E))./(Energy_vec(E+1)-Energy_vec(E));
    end
end
end