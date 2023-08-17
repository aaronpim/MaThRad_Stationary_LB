function [S] = Function_Mean_Energy_Loss(Inclusion,DX,alpha_inc,alpha_exc,p_inc,p_exc)
%##########################################################################
% Inputs
%   Inclusion   An [N,N,N,M,En] logical tensor which denotes the position
%               of the hard inclusion in the domain.
%   DX          An [M,1] double tensor which denotes the mesh length in
%               each trajectory direction.
%   alpha_inc   A positive double value that represents the coefficent of
%               the mean energy loss per unit length in the inclusion.
%   alpha_exc   A positive double value that represents the coefficent of
%               the mean energy loss per unit length in the exclusion.
%   p_inc       A positive double value that represents the power of the
%               mean energy loess per unit length in the inclusion.
%   p_exc       A positive double value that represents the power of the
%               mean energy loess per unit length in the exclusion.
% Outputs
%   S           An [N,N,N,M,En] double tensor which denotes the mean energy
%               loss per unit length.
%##########################################################################

% Define Inclusion
[~,~,~,M,En] = size(Inclusion);

% Define Energy levels
[Energy_vec] = Function_Energy_vec(En);

S_inc = @(i) max(alpha_inc*(1./Energy_vec(i))*(log(p_inc*Energy_vec(i)./(1-Energy_vec(i)))-Energy_vec(i)),0);
S_exc = @(i) max(alpha_exc*(1./Energy_vec(i))*(log(p_exc*Energy_vec(i)./(1-Energy_vec(i)))-Energy_vec(i)),0);

%S_inc = @(i) alpha_inc*(Energy_vec(i))^p_inc;
%S_exc = @(i) alpha_exc*(Energy_vec(i))^p_exc;
S = zeros(size(Inclusion));
for E = 1:En
    for i = 1:M
        S(:,:,:,i,E) = S_inc(E)*Inclusion(:,:,:,i,E)./DX(i) ...
                        + S_exc(E)*(~Inclusion(:,:,:,i,E))./DX(i);
    end
end
end