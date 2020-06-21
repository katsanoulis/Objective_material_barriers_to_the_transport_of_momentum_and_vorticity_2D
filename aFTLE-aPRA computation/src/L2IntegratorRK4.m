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
function [xt,yt] = L2IntegratorRK4(x0,y0,tspan,L2x_i,L2y_i)
x0 = x0(:); y0 = y0(:);
%% Computing the final positions of the Lagrangian particles:
    
    xt = zeros(2,length(x0));
    yt = zeros(2,length(y0));

    xt(1,:) = x0;
    yt(1,:) = y0;

    for s=1:numel(tspan)-1

        ds = tspan(s+1) - tspan(s);

        [UK1,VK1] = ODEfun(x0,y0,L2x_i,L2y_i);
        xx = x0 + 0.5 * ds * UK1;
        yy = y0 + 0.5 * ds * VK1;

        [UK2,VK2] = ODEfun(xx,yy,L2x_i,L2y_i);
        xx = x0 + 0.5 * ds * UK2;
        yy = y0 + 0.5 * ds * VK2;

        [UK3,VK3] = ODEfun(xx,yy,L2x_i,L2y_i);
        xx = x0 + ds * UK3;
        yy = y0 + ds * VK3;

        [UK4,VK4] = ODEfun(xx,yy,L2x_i,L2y_i);
        
        %increment in trajectories (displacement of grid)
        deltax = ds / 6 * (UK1 + 2 * UK2 + 2 * UK3 + UK4);
        deltay = ds / 6 * (VK1 + 2 * VK2 + 2 * VK3 + VK4);

        %update particle positions
        xx = x0 + deltax;
        yy = y0 + deltay;

        x0 = xx;
        y0 = yy;


    end
    xt(2,:) = x0;
    yt(2,:) = y0;
end

function [U1,V1] = ODEfun(y1,y2,L2x_i,L2y_i)

    U1 = L2x_i( wrapTo2Pi(y1),wrapTo2Pi(y2) );
    V1 = L2y_i( wrapTo2Pi(y1),wrapTo2Pi(y2) );

end
    