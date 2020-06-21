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

%% Define the boundaries of the domain to analyze
% resolution of points in x-y direction
mx = 350; my = 350;
xi_span = linspace(2.8,4.9,mx); yi_span = linspace(1,3,my);

%% Define parameters of the grid over which the barrier fields are defined
xmax = 6.28; ymax = 6.28;
numX = 1100; numY = 1100;

%% Compute largest rate-of-strain eigenvalue
[lambda_S2,xspan_grid,yspan_grid] = LargestRateOfStrainEigenvalue();

%% barrier field integration parameters for the instantaneous momentum barrier field
% number of dummy time steps used for integration with RK4
tSteps = 101;
% final dummy time
sFinal = 0.05;

%% Compute aFTLE and aPRA for the instantaneous momentum barrier field
[FTLEm,PRAm] = ComputeInstantaneousActiveMomentumFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY);

%% barrier field integration parameters for the instantaneous vorticity barrier field
% number of dummy time steps used for integration with RK4
tSteps = 101;
% final dummy time
sFinal = 0.15;

%% Compute aFTLE and aPRA for the instantaneous vorticity barrier field
[FTLEv,PRAv] = ComputeInstantaneousActiveVorticityFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY);

%% Comparison of FTLE and aFTLE

figure();
fontsizeaxlab = 28;
AxthicksFnt = 28;
set(gcf,'color','w');
set(gcf, 'Position', [0, 0, 1645, 445])

for fi = 1:3
if fi == 1
    subplot(1,3,fi);
    pcolor(xspan_grid,yspan_grid,lambda_S2');
    shading interp
    title('$$(a)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 max(lambda_S2(:))])
elseif fi == 2
    subplot(1,3,fi);
    pcolor(xi_span,yi_span,FTLEm');
    shading interp
    title('$$(b)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 100])
else
    subplot(1,3,fi);
    pcolor(xi_span,yi_span,FTLEv');
    shading interp
    title('$$(c)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 20])
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
    set(get(hhF,'xlabel'),'string','FTLE$$_{0}^{0}\left(\mathbf{x}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
elseif fi == 2
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aFTLE$$_{0,0}^{0.05}\left(\mathbf{x};\rho\mathbf{u}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
else
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aFTLE$$_{0,0}^{0.15}\left(\mathbf{x};\mathbf{\omega}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
end
box on
end

%% Comparison of aPRAs

hf = figure();
fontsizeaxlab = 28;
AxthicksFnt = 28;
set(gcf,'color','w');
set(gcf, 'Position', [0, 0, 1245, 445])
for fi = 1:2
if fi == 1
    subplot(1,2,fi);
    pcolor(xi_span,yi_span,PRAm');
    shading interp
    title('$$(a)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
    caxis([0 max(PRAm(:))])
else
    subplot(1,2,fi);
    pcolor(xi_span,yi_span,PRAv');
    shading interp
    title('$$(b)$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
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
    set(get(hhF,'xlabel'),'string','aPRA$$_{0,0}^{0.05}\left(\mathbf{x};\rho\mathbf{u}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
else
    hhF=colorbar(gca);
    hhF.Location='eastOutside';
    set(get(hhF,'xlabel'),'string','aPRA$$_{0,0}^{0.15}\left(\mathbf{x};\mathbf{\omega}\right)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
    set(hhF,'TickLabelInterpreter', 'latex');
end

box on
end