function plot_LB_inflow(u,logyn)
%##########################################################################
% Inputs
%   u           An [N,N,N,26,En] double tensor which represents the
%               distribution of particles in the domain.
%   logyn       A logical variable which determines if the scale is
%               logarithmic or not.
%##########################################################################

% Default values
if nargin == 1
    logyn = true;
end

% Define the position of the boundary.
[N,M,En] = size(u,[1,4,5]);
[~,~,Gamma_mnus] = Boundary_Gen(N,En);

% Define a temporary tensor which we assign to the boundary
V = zeros(size(u));
V(Gamma_mnus) = u(Gamma_mnus);
V = sum(V,[4,5]);
if logyn
    V = log10(sum(V,[4,5]));    
end

% Define 3D co-ordinates
[X,Y,Z] = meshgrid(linspace(-1,1,N),linspace(-1,1,N),linspace(-1,1,N));

% Plot along each of the faces
figure; hold on
h = slice(X,Y,Z,V,[-1],[],[]); h.FaceAlpha = 0.5; % x =-1
h = slice(X,Y,Z,V,[ 1],[],[]); h.FaceAlpha = 0.5; % x = 1
h = slice(X,Y,Z,V,[],[-1],[]); h.FaceAlpha = 0.5; % y =-1
h = slice(X,Y,Z,V,[],[ 1],[]); h.FaceAlpha = 0.5; % y = 1
h = slice(X,Y,Z,V,[],[],[-1]); h.FaceAlpha = 0.5; % z =-1
h = slice(X,Y,Z,V,[],[],[ 1]); h.FaceAlpha = 0.5; % z = 1

% Change visual settings
shading interp; axis equal; cbh = colorbar;

% If log scale then redo the colorbar ticks.
if logyn
    cbh.TickLabels = 10.^(cbh.Ticks);
end

% Plot outline of domain
plot3([-1,1,1,-1,-1],[-1,-1,1,1,-1],[-1,-1,-1,-1,-1],'-k')
plot3([-1,1,1,-1,-1],[-1,-1,1,1,-1],[ 1, 1, 1, 1, 1],'-k')
plot3([-1,-1],[-1,-1],[-1, 1],'-k'); plot3([1,1],[-1,-1],[-1, 1],'-k');
plot3([-1,-1],[ 1, 1],[-1, 1],'-k'); plot3([1,1],[ 1, 1],[-1, 1],'-k');

% Define dense inclusion
A = Inclusion_defn(N,M,size(u,5));
A = A(:,:,:,1,1);

% Plot outline of dense inclusion
s = isosurface(X,Y,Z,A,0.99); p = patch(s); isonormals(X,Y,Z,A,p);

% Change visual settings
set(p,'FaceColor',[1 0.5 0.5],'FaceAlpha',0.0,'EdgeColor',[0 0 0]);

% Define the target
[~,Target] = Target_Function(u,0); B = Target(:,:,:,1);

% Plot outline of the target
s = isosurface(X,Y,Z,B,0.99); p = patch(s); isonormals(X,Y,Z,B,p);

% Change visual settings
view(210,45)
set(p,'FaceColor',[1 0.5 0.5],'FaceAlpha',0.0,'EdgeColor',[1 0 0]); 
hold off
end