% This runs the forward solver for given incident energy beam, and plots
% the associated depth dose curves for varying width.

% Define coefficients
N = 50;
En= 50;

% Define trajectory vectors and boundary
[s_vec,~,Gamma_mnus] = Boundary_Gen(N,En);

% Define the dense inclusion
M = length(s_vec);
[Inclusion] = Function_Inclusion_defn(N,M,En);

% Define Boundary conditions
BC = zeros(N,N,N,length(s_vec),En);
if ~logical(mod(N,2))
    mid = N/2:N/2+1;
    BC(1,mid,mid,sum(abs(s_vec),2)==1, round(En/2)) = 1;
else
    mid = (N-1)/2:(N+3)/2;
    BC(1,mid,mid,sum(abs(s_vec),2)==1, round(En/2)) = 1;
end
Domain_width_vec = 10.^(linspace(-2,2,50)); % Measured in metres^-1

% Define absorpsion & scattering tensor
[sigma_a] = Function_sigma_a_dist(Inclusion);
[sigma_s] = Function_sigma_s_dist(Inclusion, s_vec);
      
for i = 1:length(Domain_width_vec)
    Domain_width = Domain_width_vec(i);
    [S] = Function_Mean_Energy_Loss(Inclusion,(2/(N-1))*sqrt(sum(s_vec.^2,2)),Domain_width*1.2658e-5, Domain_width*1.2658e-5, 13626.64, 13626.64);
    % Run forward
    [u,uerror] = LB_primal(Gamma_mnus, S, sigma_s, sigma_a, s_vec, BC, 100, 1e-7);
    A = Function_Dose_Calculation(u); A = A(:,mid,mid); A = sum(A,[2,3]);
    save(['Test_forward_width_vary_energy_25_',num2str(i),'.mat'],"A","u","uerror")
end
for i = 1:50
    load(['Test_forward_width_vary_energy_25_',num2str(i),'.mat'],"A")
    B = max(A,(1.0e-8)*ones(size(A)));
    figure; plot(linspace(-1,1,50),B)
end