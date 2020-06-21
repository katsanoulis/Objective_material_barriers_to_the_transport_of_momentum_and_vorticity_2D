function [FTLE,x0,y0] = ClassicFTLE(xt,yt,tc)
%% FTLE computation for velocity field

x0 = squeeze(xt(1,:,:,5));
y0 = squeeze(yt(1,:,:,5));

rhox = xt(1,1,1,1) - xt(1,1,1,5);
rhoy = rhox;

F11 = squeeze((xt(tc,:,:,1)-xt(tc,:,:,3)))/(2*rhox);
F12 = squeeze((xt(tc,:,:,2)-xt(tc,:,:,4)))/(2*rhox);
F21 = squeeze((yt(tc,:,:,1)-yt(tc,:,:,3)))/(2*rhoy);
F22 = squeeze((yt(tc,:,:,2)-yt(tc,:,:,4)))/(2*rhoy);

C11 = F11.^2+F21.^2;
C12 = F11.*F12+F22.*F21;
C22 = F22.^2+F12.^2;

trC  = C11+C22;
detC = C11.*C22-C12.^2;

lambda_2 = 0.5*trC+sqrt((0.5*trC).^2-detC);

FTLE = 1/2/(25)*log(lambda_2);

end