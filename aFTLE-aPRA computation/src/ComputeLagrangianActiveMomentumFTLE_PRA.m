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
function [FTLEm,PRAm] = ComputeLagrangianActiveMomentumFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY,xt,yt,tc)
%% Load grid vectors
xspan_over = linspace(0,xmax,numX);
yspan_over = linspace(0,ymax,numY);

[xgrid_over,ygrid_over] = ndgrid(xspan_over(1:end-1),yspan_over(1:end-1));

% Grid over which velocity is defined
load('../../Data/turb_w_0000.mat','x')
xspan_grid = x(1,:);
xspan_grid = [xspan_grid 2*pi];
yspan_grid = xspan_grid;

[m,n] = size(xgrid_over);

mx = length(xi_span); my = length(yi_span);

%% Compute time-averaged Lagrangian vorticity
tspan = 0:tc;
deltat = tc - tspan(1);
dt = (deltat)/(length(tspan) - 1);

str1 = '../../Data/turb_w_';
str2 = pad(int2str(0),4,'left','0');
str = strcat(str1,str2);
load(str);

Lomega = zeros(m,n);

w = [w w(:,1)]; w=[w; w(1,:)]';

omega = griddedInterpolant({xspan_grid,yspan_grid},w,'spline','none');

omegaInst = omega(wrapTo2Pi(xgrid_over),wrapTo2Pi(ygrid_over));

Lomega = Lomega + omegaInst*dt/deltat/2;

for num = 2:length(tspan)-1

    str1 = '../../Data/turb_w_';
    u_int = 5*num;
    str2 = pad(int2str(u_int),4,'left','0');
    str = strcat(str1,str2);
    load(str);
    
    w = [w w(:,1)]; w=[w; w(1,:)]';
    
    omega = griddedInterpolant({xspan_grid,yspan_grid},w,'spline','none');
    
    x0 = squeeze(xt(num,:,:,5));
    y0 = squeeze(yt(num,:,:,5));

    omegaInst = omega(wrapTo2Pi(x0),wrapTo2Pi(y0));

    Lomega = Lomega + omegaInst*dt/deltat;
    
end

num = length(tspan);
str1 = '../../Data/turb_w_';
u_int = 5*num;
str2 = pad(int2str(u_int),4,'left','0');
str = strcat(str1,str2);
load(str);

w = [w w(:,1)]; w=[w; w(1,:)]';
    
omega = griddedInterpolant({xspan_grid,yspan_grid},w,'spline','none');

x0 = squeeze(xt(num,:,:,5));
y0 = squeeze(yt(num,:,:,5));

omegaInst = omega(wrapTo2Pi(x0),wrapTo2Pi(y0));

Lomega = Lomega + omegaInst*dt/deltat/2;

Lomega = [Lomega Lomega(:,1)];
Lomega = [Lomega; Lomega(1,:)];

%% Compute interpolants for the Lagrangian momentum barrier field

dx = abs( xspan_over(2)-xspan_over(1) );
dy = abs( yspan_over(2)-yspan_over(1) );
[by,bx] = gradient(Lomega',dx,dy);
bx = -bx;

bx_interp_m = griddedInterpolant({[xspan_over(1:end-1) 2*pi],[yspan_over(1:end-1) 2*pi]},bx','spline','none');
by_interp_m = griddedInterpolant({[xspan_over(1:end-1) 2*pi],[yspan_over(1:end-1) 2*pi]},by','spline','none');

%% Integration of momentum barrier field

[xi,yi] = ndgrid(xi_span,yi_span);
Nrad = 4;

% Define auxiliary grid for the computation of the active flowmap gradients
rhox = (xspan_over(2) - xspan_over(1))*0.01;
rhoy = rhox;

xit=zeros(mx,my,Nrad);
yit=zeros(mx,my,Nrad);
for k=1:Nrad
    xit(:,:,k) = xi + rhox*cos( (k-1)*pi/2 );
    yit(:,:,k) = yi + rhoy*sin( (k-1)*pi/2 );
end

sSpan = linspace(0,sFinal,tSteps);

[xit,yit] = L2IntegratorRK4(xit,yit,sSpan,bx_interp_m,by_interp_m);

%% Compute the active momentum FTLE
xit = reshape(xit, 2, mx, my, Nrad);
yit = reshape(yit, 2, mx, my, Nrad);

tcb = 2;

F11 = squeeze((xit(tcb,:,:,1)-xit(tcb,:,:,3)))/(2*rhox);
F12 = squeeze((xit(tcb,:,:,2)-xit(tcb,:,:,4)))/(2*rhox);
F21 = squeeze((yit(tcb,:,:,1)-yit(tcb,:,:,3)))/(2*rhoy);
F22 = squeeze((yit(tcb,:,:,2)-yit(tcb,:,:,4)))/(2*rhoy);

C11 = F11.^2+F21.^2;
C12 = F11.*F12+F22.*F21;
C22 = F22.^2+F12.^2;

trC  = C11+C22;
detC = C11.*C22-C12.^2;

lambda_2 = 0.5*trC+sqrt((0.5*trC).^2-detC);

FTLEm = 1/2/(sSpan(end) - sSpan(1))*log(lambda_2);
%FTLEm = 1/2/(1 - sSpan(1))*log(lambda_2);

%% Compute the active momentum PRA

PRAm = zeros(mx,my);
for i = 1:mx
	for j = 1:my
		DF = [F11(i,j) F12(i,j)
		      F21(i,j) F22(i,j)];
		[U,S,V] = svd(DF);
		
		ksi1 = squeeze(V(:,1));
		%ksi2 = squeeze(V(:,2));

		eta1 = squeeze(U(:,1));
		eta2 = squeeze(U(:,2));

        si = ksi1(1)*eta2(1) + ksi1(2)*eta2(2);
        co = ksi1(1)*eta1(1) + ksi1(2)*eta1(2);
        PRAm(i,j) = (1-sign(si))*pi + sign(si)*acos(co);

	end
end

PRAm = real(PRAm);
end

