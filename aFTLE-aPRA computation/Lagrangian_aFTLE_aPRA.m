%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%              ____________________   ___           %%%%%%%%%%%%%
%%%%%%%%%%%             /  ________   ___   /__/  /           %%%%%%%%%%%%%
%%%%%%%%%%%            /  _____/  /  /  /  ___   /            %%%%%%%%%%%%%
%%%%%%%%%%%           /_______/  /__/  /__/  /__/             %%%%%%%%%%%%%
%%%%%%%%%%%    Swiss Federal Institute of Technology Zurich   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Stergios Katsanoulis  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Email:  katsanos@ethz.ch      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Date:   18/06/2020            %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1]: G. Haller, S. Katsanoulis, M. Holzner, B. Frohnapfel & D. Gatti,
%      Objective material barriers to the transport of momentum and vorticity. submitted (2020)
%% Add necessary folders to MATLAB path
addpath('src')

%% Load advected positions of initial grid
load('../../particles1100.mat')

% Specify final time to analyze: [t0,t1] = [0,tc]
tc = 25;

%% Define the boundaries of the domain to analyze
% resolution of points in x-y direction
mx = 550; my = 550;
xi_span = linspace(2.8,4.9,mx); yi_span = linspace(1,3,my);

%% Define parameters of the grid over which the barrier fields are defined
xmax = 6.28; ymax = 6.28;
numX = 1100; numY = 1100;

%% Compute Classic FTLE
[FTLE,x0,y0] = ClassicFTLE(xt,yt,tc);

%% Compute Classic PRA
PRA = ClassicPRA(xt,yt,tc);

%% barrier field integration parameters for the Lagrangian momentum barrier field
% number of dummy time steps used for integration with RK4
tSteps = 701;
% final dummy time
sFinal = 0.35;

%% Compute aFTLE and aPRA for the Lagrangian momentum barrier field
[FTLEm,PRAm] = ComputeLagrangianActiveMomentumFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY,xt,yt,tc);

%% barrier field integration parameters for the Lagrangian vorticity barrier field
% number of dummy time steps used for integration with RK4
tSteps = 1001;
% final dummy time
sFinal = 0.05;

%% Compute aFTLE and aPRA for the Lagrangian vorticity barrier field
[FTLEv,PRAv] = ComputeLagrangianActiveVorticityFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY,xt,yt,tc);

%% Comparison of FTLE and aFTLE
figure();
fontsizeaxlab = 28;
AxthicksFnt = 28;
set(gcf,'color','w');
set(gcf, 'Position', [0, 0, 1645, 445])
for fi = 1:3

if fi == 1
    subplot(1,3,fi);
    pcolor(x0,y0,FTLE);
    shading interp
    title('$$(a)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 max(FTLE(:))])
elseif fi == 2
    subplot(1,3,fi);
    pcolor(xi_span,yi_span,FTLEm');
    shading interp
    title('$$(b)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 25])
else
    subplot(1,3,fi);
    pcolor(xi_span,yi_span,FTLEv');
    shading interp
    title('$$(c)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 100])
end

set(gca,'YDir','normal')
axis equal
xlabel('$$x$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
if fi==1
ylabel('$$y$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
end
axis([2.8 4.9 1 3]);
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',AxthicksFnt,'fontWeight','normal');

if fi == 1
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','FTLE$$_{0}^{25}\left(\mathbf{x}_{0}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
elseif fi == 2
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aFTLE$$_{0,25}^{0.35}\left(\mathbf{x}_{0};\rho\mathbf{u}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
else
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aFTLE$$_{0,25}^{0.05}\left(\mathbf{x}_{0};\mathbf{\omega}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
end
box on
end

%% Comparison of PRA and aPRA
hf = figure();
fontsizeaxlab = 28;
AxthicksFnt = 28;
set(gcf,'color','w');
set(gcf, 'Position', [0, 0, 1645, 445])
for fi = 1:3

if fi == 1
    subplot(1,3,fi);
    pcolor(x0,y0,PRA);
    shading interp
    title('$$(a)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 max(PRA(:))])
elseif fi == 2
    subplot(1,3,fi);
    pcolor(xi_span,yi_span,PRAm');
    shading interp
    title('$$(b)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 max(PRAm(:))])
else
    subplot(1,3,fi);
    pcolor(xi_span,yi_span,PRAv');
    shading interp
    title('$$(c)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 max(PRAv(:))])
end

set(gca,'YDir','normal')
axis equal
xlabel('$$x$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
if fi==1
ylabel('$$y$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
end
axis([2.8 4.9 1 3]);
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',AxthicksFnt,'fontWeight','normal');

if fi == 1
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','PRA$$_{0}^{25}\left(\mathbf{x}_{0}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
elseif fi == 2
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aPRA$$_{0,25}^{0.35}\left(\mathbf{x}_{0};\rho\mathbf{u}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
else
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aPRA$$_{0,25}^{0.05}\left(\mathbf{x}_{0};\mathbf{\omega}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
end
box on
end