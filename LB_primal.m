function [u,uerror] = LB_primal(Gamma_mnus, S, sigma_s, sigma_a, s_vec, BC, n_iter,min_error)
%% This code computes a numerical solution to the Boltzmann-CSD equation.
%##########################################################################
% Inputs
%   sigma_s     An [N,N,N,26,26,En] double tensor that represents the 
%               probability that a particle at a given position and energy 
%               scatters with a particlar angle.
%   sigma_a     An [N,N,N,26,En] double tensor which denotes the absorpsion
%               coefficient at each point, each trajectory and each
%               energy.
%   BC          An [N,N,N,26,En] double tensor which represents the inflow
%               boundary data
%   n_iter      An integer that is the maximum number of iterations before
%               termination
%   min_error   A positive double that represents the smallest value the
%               relative error can take  before termination.
% Outputs
%   u           An [N,N,N,26,En] double tensor which represents the
%               distribution of particles in the domain.
%   uerror      An [n_iter,1] vector which stores the relative error
%               between each iteration.
%##########################################################################
%% Preallocation functions
u = zeros(size(Gamma_mnus));
L2 = @(V)   sqrt(sum(V.^2,"all"));   % L2 norm
    
% Apply the inflow BC's and normalise
u(Gamma_mnus) = BC(Gamma_mnus); 
u = u./sum(Function_Dose_Calculation(u),"all");

% Termination criteria
count = 0; uerror = zeros(n_iter+1,1); uerror(1) = inf;

while count < n_iter && uerror(count) > min_error
    % Update termination criteria
    count = count +1;   uold = u;

    % Calculate particle scattering
    Scatter_Ten = Scattering(u,sigma_s,s_vec)-sigma_a.*u+Energy_loss(u,S);
    % Scattering = int( sigma_s u ) ds'
    % Energy_loss= d/dE( S u )

    % Calculate particle streaming
    [u] = Streaming(u, Scatter_Ten, s_vec);
    
    % Enforce the boundary conditions
    u(Gamma_mnus) = BC(Gamma_mnus);  

    % Remove any points of negative flux.
    u(u<0)=0;
    
    % Remove any particles of zero energy.
    u(:,:,:,:,1)=0;

    % Normalise
    u = u./sum(Function_Dose_Calculation(u),"all");
    
    % Increment error
    uerror(count+1) = L2(u-uold)./L2(u);

end
end
