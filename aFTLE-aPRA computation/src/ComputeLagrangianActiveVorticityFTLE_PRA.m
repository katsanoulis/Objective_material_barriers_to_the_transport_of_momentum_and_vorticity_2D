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
function [FTLEv,PRAv] = ComputeLagrangianActiveVorticityFTLE_PRA(xi_span,yi_span,tSteps,sFinal,xmax,ymax,numX,numY,xt,yt,tc)
%% Load grid vectors
xspan_over = linspace(0,xmax,numX);
yspan_over = linspace(0,ymax,numY);

[xgrid_over,ygrid_over] = ndgrid(xspan_over(1:end-1),yspan_over(1:end-1));

% Grid over which velocity is defined
load('../../Data/turb_w_0000.mat','x')
xspan_grid = x(1,:);
xspan_grid = [xspan_grid 2*pi];
yspan_grid = xspan_grid;

mx = length(xi_span); my = length(yi_span);

%% Compute Lagrangian vorticity-change function

t0 = 0;
t1 = tc;

load('../../Data/turb_w_0000.mat','w')
w = [w w(:,1)];
w = [w; w(1,:)]';

omega0 = griddedInterpolant({xspan_grid,yspan_grid},w,'spline','none');

str1 = '../../Data/turb_w_';
u_int = 5*tc+1;
str2 = pad(int2str(u_int),4,'left','0');
str = strcat(str1,str2);
load(str);
w = [w w(:,1)];
w = [w; w(1,:)]';

omega1 = griddedInterpolant({xspan_grid,yspan_grid},w,'spline','none');

w0 = omega0(wrapTo2Pi(xgrid_over),wrapTo2Pi(ygrid_over));

x1 = squeeze(xt(tc+1,:,:,5)); y1 = squeeze(yt(tc+1,:,:,5));
w1 = omega1(wrapTo2Pi(x1),wrapTo2Pi(y1));

%vIncrement = (w1-w0)/(t1-t0);
vIncrement = w1-w0;

vIncrement = [vIncrement vIncrement(:,1)];
vIncrement = [vIncrement; vIncrement(1,:)];

%% Compute interpolants for the Lagrangian vorticity barrier field

dx = abs( xspan_grid(2)-xspan_grid(1) );
dy = abs( yspan_grid(2)-yspan_grid(1) );
[by,bx] = gradient(vIncrement',dx,dy);
bx = -bx;

bx_interp_v = griddedInterpolant({[xspan_over(1:end-1) 2*pi],[yspan_over(1:end-1) 2*pi]},bx','spline','none');
by_interp_v = griddedInterpolant({[xspan_over(1:end-1) 2*pi],[yspan_over(1:end-1) 2*pi]},by','spline','none');

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

