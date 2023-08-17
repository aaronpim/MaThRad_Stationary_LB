function [scatter_ten] = Scattering(u,sigma_s,s_vec)
%##########################################################################
% Inputs
%   sigma_s     An [N,N,N,26,26,En] double tensor that represents the 
%               probability that a particle at a given position and energy 
%               scatters with a particlar angle.
%   u           An [N,N,N,26,En] double tensor that represents the particle
%               probability distribution in terms of position, trajectory 
%               and energy.
% Outputs
%   scatter_ten An [N,N,N,26,En] double tensor that represents the integral
%               of sigma_s and u with respect to trajectory.
%##########################################################################

% Pre-allocate scatter_ten
scatter_ten = zeros(size(u));   
M = size(u,4);
[~,assign_matrix] = dot_matrix(s_vec);

% Compute the tensor multiplication
for i = 1:M %Exit trajectory
    for j = 1:M %Entry trajectory  
        scatter_ten(:,:,:,i,:) = scatter_ten(:,:,:,i,:)+ sigma_s(:,:,:,assign_matrix(i,j),:).*u(:,:,:,j,:);
    end
end

% Normalise
scatter_ten = scatter_ten./sum(scatter_ten,"all");
end