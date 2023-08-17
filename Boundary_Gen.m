function [s_vec,Gamma_plus,Gamma_mnus] = Boundary_Gen(N,En)
%##########################################################################
% Inputs
%   N           The number of spatial nodes along each axis.
%   En          The number of energy levels.
% Outputs
%   s_vec       A 26x3 matrix denotes the possible trajectories a partcle 
%               can take.
%   Gamma_plus  An [N,N,N,26,En] logical tensor which denotes the position
%               of the outflow boundary.
%   Gamma_mnus  An [N,N,N,26,En] logical tensor which denotes the position
%               of the inflow boundary.
%##########################################################################
%%  Define trajectory vector
s_vec = combvec([-1,0,1],[-1,0,1],[-1,0,1])'; 
s_vec = s_vec(sum(s_vec.^2,2)~=0,:);
M = length(s_vec);

%% Preallocate the inflow & outflow boundaries
Gamma_plus = false(N,N,N,M,En);    
Gamma_mnus = false(N,N,N,M,En);
[X,Y,Z] = meshgrid(1:N,1:N,1:N);
X = permute(X,[2,1,3]);
Y = permute(Y,[2,1,3]);
a = 2:N-1;
%% Define Gamma_mnus & Gamma_plus
for s = 1:length(s_vec)
    for E = 1:En
        % Project (X,Y,Z) into (X,Y,Z)+s(i,:).
        
        %   If the transformed cordinates are in the interior of the domain
        %   then they are an element of Gamma_mnus.
        Gamma_mnus(:,:,:,s,E) = ismember(X+s_vec(s,1),X(a,a,a)) & ismember(Y+s_vec(s,2),Y(a,a,a)) & ismember(Z+s_vec(s,3),Z(a,a,a));
        
        %   If the transformed coordinates are not contained in the domain
        %   then they are an element of Gamma_plus
        Gamma_plus(:,:,:,s,E) = ~ismember(X+s_vec(s,1),X) | ~ismember(Y+s_vec(s,2),Y) | ~ismember(Z+s_vec(s,3),Z);
    end
end
% Delete any elements of Gamma_plus & Gamma_mnus that are not on the
% boundary.
Gamma_mnus(a,a,a,:,:) = false;
Gamma_plus(a,a,a,:,:) = false;
end