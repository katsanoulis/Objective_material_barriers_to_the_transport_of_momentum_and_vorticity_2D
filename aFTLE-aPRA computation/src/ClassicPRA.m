function PRA = ClassicPRA(xt,yt,tc)
x0 = squeeze(xt(1,:,:,5));
y0 = squeeze(yt(1,:,:,5));

rhox = xt(1,1,1,1) - xt(1,1,1,5);
rhoy = rhox;

F11 = squeeze((xt(tc,:,:,1)-xt(tc,:,:,3)))/(2*rhox);
F12 = squeeze((xt(tc,:,:,2)-xt(tc,:,:,4)))/(2*rhox);
F21 = squeeze((yt(tc,:,:,1)-yt(tc,:,:,3)))/(2*rhoy);
F22 = squeeze((yt(tc,:,:,2)-yt(tc,:,:,4)))/(2*rhoy);

PRA = 0*x0;
for i = 1:size(x0,1)
	for j = 1:size(x0,2)
		DF = [F11(i,j) F12(i,j)
		      F21(i,j) F22(i,j)];
		[U,S,V] = svd(DF);
		
		ksi1 = squeeze(V(:,1));
		ksi2 = squeeze(V(:,2));

		eta1 = squeeze(U(:,1));
		eta2 = squeeze(U(:,2));

        si = ksi1(1)*eta2(1) + ksi1(2)*eta2(2);
        co = ksi1(1)*eta1(1) + ksi1(2)*eta1(2);
        PRA(i,j) = (1-sign(si))*pi + sign(si)*acos(co);

	end
end

end