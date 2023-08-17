function [sigma_s] = Function_sigma_s_dist(Inclusion, s_vec)
%##########################################################################
% Inputs
%   Inclusion   An [N,N,N,M,En] logical tensor which denotes the position
%               of the dense inclusion.
%   s_vec       An Mx3 tensor that denotes the possible trajectories a
%               particle can take.
% Outputs
%   sigma_s     An [N,N,N,M,M,En] double tensor that represents the 
%               probability that a particle at a given position and energy 
%               scatters with a particlar angle.
%##########################################################################

    [N,En] = size(Inclusion,[1,5]);

% Define dot product space
    [vals,~] = dot_matrix(s_vec);

% Pre-allocate the sigma_s vector
    sigma_s = zeros(N,N,N,length(vals),En);

% Define Rutherford Scattering kernal
    R = @(d,a) 1./((1+a-d));

% The parameter a>0 controls how highly forwardly peaked the scattering is.
 
% Define how the parameter a varies with respect to energy in the inclusion
    a_inc = @(E) E.^(-0.2);

% Define how the parameter a varies with respect to energy in the exclusion
    a_exc = @(E) 0.0001*E.^(-1);

% Define the energy vec
    [Energy_vec] = Function_Energy_vec(En);

% Define sigma_s
for i = 1:length(vals)
    for j = 1:En    % Particle energy
        sigma_s(:,:,:,i,j) = Inclusion(:,:,:,i,j)*R(vals(i),a_inc(Energy_vec(j)))...
            + (~Inclusion(:,:,:,i,j))*R(vals(i),a_exc(Energy_vec(j)));
    end
end

% Warning sigma_s is NOT normalised, the normalisation step will be applied
% during the scattering operation.
end