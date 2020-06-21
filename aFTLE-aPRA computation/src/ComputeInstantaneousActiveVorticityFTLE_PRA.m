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
function [FTLEv,PRAv] = ComputeInstantaneousActiveVorticityFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY)
%% Load grid vectors
xspan_over = linspace(0,xmax,numX);
yspan_over = linspace(0,ymax,numY);

% Grid over which velocity is defined
load('../../Data/turb_w_0000.mat','x')
xspan_grid = x(1,:);
xspan_grid = [xspan_grid 2*pi];
yspan_grid = xspan_grid;

[x0_grid,y0_grid] = ndgrid(xspan_grid,yspan_grid);

mx = length(xi_span); my = length(yi_span);

%% Compute instantaneous acceleration
tspan = [0 0.1 0.2];
    
[m,n] = size(x0_grid);
xt=x0_grid;
yt=y0_grid;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

tic
[xt,yt] = Integrator_velocity(xt(:,:),yt(:,:),tspan,options);

disp('Particle advection is over.')
toc

xt = reshape(xt, length(tspan), m, n);
yt = reshape(yt, length(tspan), m, n);

Ax = squeeze((xt(3,:,:)-2*xt(2,:,:)+xt(1,:,:))/tspan(2)^2);
Ay = squeeze((yt(3,:,:)-2*yt(2,:,:)+yt(1,:,:))/tspan(2)^2);

%% Compute curl of instantaneous acceleration
Ax = permute(Ax,[2 1 3]);
Ay = permute(Ay,[2 1 3]);

[x0_curl,y0_curl] = meshgrid(xspan_grid,yspan_grid);

Cz = curl(x0_curl,y0_curl,Ax,Ay);

%% Compute interpolants for the instantaneous vorticity barrier field

dx = abs( xspan_grid(2)-xspan_grid(1) );
dy = abs( yspan_grid(2)-yspan_grid(1) );
[by,bx] = gradient(Cz,dx,dy);
bx = -bx;

bx_interp_v = griddedInterpolant({xspan_grid,yspan_grid},bx','spline','none');
by_interp_v = griddedInterpolant({xspan_grid,yspan_grid},by','spline','none');

%% Integration of vorticity barrier field
[xi,yi] = ndgrid(xi_span,yi_span);
Nrad = 4;

% Define auxiliary grid for the computation of the active flowmap gradients
rhox = (xspan_over(2) - xspan_over(1))*0.01;
rhoy = rhox;

xvt=zeros(mx,my,Nrad);
yvt=zeros(mx,my,Nrad);
for k=1:Nrad
    xvt(:,:,k) = xi + rhox*cos( (k-1)*pi/2 );
    yvt(:,:,k) = yi + rhoy*sin( (k-1)*pi/2 );
end

sSpan = linspace(0,sFinal,tSteps);

[xvt,yvt] = L2IntegratorRK4(xvt,yvt,sSpan,bx_interp_v,by_interp_v);

%% Compute the active momentum FTLE
xvt = reshape(xvt, 2, mx, my, Nrad);
yvt = reshape(yvt, 2, mx, my, Nrad);

tc = 2;

F11 = squeeze((xvt(tc,:,:,1)-xvt(tc,:,:,3)))/(2*rhox);
F12 = squeeze((xvt(tc,:,:,2)-xvt(tc,:,:,4)))/(2*rhox);
F21 = squeeze((yvt(tc,:,:,1)-yvt(tc,:,:,3)))/(2*rhoy);
F22 = squeeze((yvt(tc,:,:,2)-yvt(tc,:,:,4)))/(2*rhoy);

C11 = F11.^2+F21.^2;
C12 = F11.*F12+F22.*F21;
C22 = F22.^2+F12.^2;

trC  = C11+C22;
detC = C11.*C22-C12.^2;

lambda_2 = 0.5*trC+sqrt((0.5*trC).^2-detC);

FTLEv = 1/2/(sSpan(end) - sSpan(1))*log(lambda_2);

%% Compute the active momentum PRA

PRAv = zeros(mx,my);
for i = 1:mx
	for j = 1:my
		DF = [F11(i,j) F12(i,j)
		      F21(i,j) F22(i,j)];
		[U,S,V] = svd(DF);
		
		ksi1 = squeeze(V(:,1));

		eta1 = squeeze(U(:,1));
		eta2 = squeeze(U(:,2));

        si = ksi1(1)*eta2(1) + ksi1(2)*eta2(2);
        co = ksi1(1)*eta1(1) + ksi1(2)*eta1(2);
        PRAv(i,j) = (1-sign(si))*pi + sign(si)*acos(co);

	end
end

PRAv = real(PRAv);
end

