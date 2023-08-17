function plot_LB(u,logyn)
%##########################################################################
% Inputs
%   u           An [N,N,N,26,En] double tensor which represents the
%               distribution of particles in the domain.
%   logyn       A logical value that determines if the dose scale is
%               logarithmic or not.
%   kappa       An [N,N,N,M,En] double tensor which denotes the stopping
%               powers of the particles in a given medium.
%##########################################################################

% Default values
if nargin == 1
    logyn = false;
end

% Define [N,N,N] dimensional spatial grid
[N,M] = size(u,[1,4]);

[X,Y,Z] = meshgrid(linspace(-1,1,N),linspace(-1,1,N),linspace(-1,1,N));

% Define Dose: V
V = Function_Dose_Calculation(u);
if logyn
    V = log10(V);
end

% Plot slices of the Dose in two planes.
figure; slice(X,Y,Z,V,[0],[],[0]); view(45,30)
shading interp; axis equal;

% Define colorbar & change ticks if logarithmic 
cbh = colorbar;
if logyn
    cbh.TickLabels = 10.^(cbh.Ticks);
end

% Plot black outline of the edges of the domain.
hold on
plot3([-1,1,1,-1,-1],[-1,-1,1,1,-1],[-1,-1,-1,-1,-1],'-k')
plot3([-1,1,1,-1,-1],[-1,-1,1,1,-1],[1,1,1,1,1],     '-k')
plot3([-1,-1],[-1,-1],[-1,1],'-k'); plot3([1,1],[-1,-1],[-1,1],'-k')
plot3([-1,-1],[1,1],[-1,1],  '-k'); plot3([1,1],[1,1],[-1,1],  '-k')

% Define inclusion Area: A
A = Function_Inclusion_defn(N,M,size(u,5)); A = A(:,:,:,1,1);

% Plot black mesh outline of the inclusion.
s = isosurface(X,Y,Z,A,0.99); p = patch(s); isonormals(X,Y,Z,A,p);

% Change visual settings
set(p,'FaceColor',[1 0.5 0.5],'FaceAlpha',0.0,'EdgeColor',[0 0 0]); 

% % Define the Target Area: Target
% [~,Target] = Target_Function(u,0);
% % Plot red mesh outline of the target volume.
% B = Target(:,:,:,1);
% s = isosurface(X,Y,Z,B,0.99);
% p = patch(s);
% isonormals(X,Y,Z,B,p)
% set(p,'FaceColor',[1 0.5 0.5]); 
% set(p,'FaceAlpha',0.0);
% set(p,'EdgeColor',[1 0 0]); 
hold off
end