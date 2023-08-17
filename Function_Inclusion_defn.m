function [Inclusion] = Function_Inclusion_defn(N,M,En)
%##########################################################################
% Inputs
%   N           The number of spatial nodes along each axis.
%   En          The number of energy levels.
%   M           The number of trajectory vectors
% Outputs
%   Inclusion   An [N,N,N,M,En] logical tensor which denotes the position
%               of the hard inclusion in the domain.    
%##########################################################################
% Define spatial co-ordinates
[X,Y,Z]= meshgrid(linspace(-1,1,N),linspace(-1,1,N),linspace(-1,1,N));

% Construct a logical [N,N,N] tensor that denotes the spatial position of
% the inclusion
Inclusion = (sqrt(X.^2+Y.^2+Z.^2)<0.7).*(sqrt(X.^2+Y.^2+Z.^2)>0.3).*(Y>0);

% Repeat the tensor until it is the correct size.
Inclusion = repmat(Inclusion,[1,1,1,M,En]);
end