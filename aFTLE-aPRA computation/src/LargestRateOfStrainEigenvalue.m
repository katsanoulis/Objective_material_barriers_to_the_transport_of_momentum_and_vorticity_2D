function [lambda_S2,xspan_grid,yspan_grid] = LargestRateOfStrainEigenvalue()
%% Load grid vectors

% Grid over which velocity is defined
load('../../Data/turb_w_0000.mat','x')
xspan_grid = x(1,:);
xspan_grid = [xspan_grid 2*pi];
yspan_grid = xspan_grid;

%%
load('../../Data/turb_u_0000.mat','u1','u2')

u1 = [u1 u1(:,1)]; u1=[u1; u1(1,:)];
u2 = [u2 u2(:,1)]; u2=[u2; u2(1,:)];

dx = abs( xspan_grid(2)-xspan_grid(1) );
dy = abs( yspan_grid(2)-yspan_grid(1) );
[ux,uy] = gradient(u1,dx,dy);
[vx,vy] = gradient(u2,dx,dy);

S11 = ux';
S12 = (uy + vx)'/2;
S22 = vy';

trS  = S11+S22;
detS = S11.*S22-S12.^2;

lambda_S2 = 0.5*trS+sqrt((0.5*trS).^2-detS);

end