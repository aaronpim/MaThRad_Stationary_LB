function [u] = Streaming(u, Scatter_ten, s_vec)
%##########################################################################
% Inputs
%   u           An [N,N,N,26,En] double tensor that represents the particle
%               probability distribution in terms of position, trajectory 
%               and energy.
%   Scatter_ten An [N,N,N,26,En] double tensor that represents the change
%               in particle energy and trajectory as the result of a
%               scattering event.
%   s_vec       A 26x3 matrix denotes the possible trajectories a partcle 
%               can take.
% Outputs
%   u           A transformed input, where all particles with trajectory S
%               move one unit in the S direction.
%##########################################################################

[N,~,~,M,En] = size(u);    s_max = max(s_vec,[],"all");

% Define Meshwidth
DX = (2/(N-1))*sqrt(sum(s_vec.^2,2));

% Define a temporary tensor
% The tensor is such that for the non-zero values of utemp, utemp(x+S) is
% still in the tensor utemp.
utemp = zeros(N+2*s_max,N+2*s_max,N+2*s_max,M,En);
a = (1+s_max):(N+s_max);
utemp(a,a,a,:,:) = u;

% Define upwind finite difference derivative.
for i = 1:M
    for E = 1:En
        utemp(a+s_vec(i,1),a+s_vec(i,2),a+s_vec(i,3),i,E) = utemp(a,a,a,i,E) + DX(i)*Scatter_ten(:,:,:,i,E);
    end
end
u = utemp(a,a,a,:,:);
end