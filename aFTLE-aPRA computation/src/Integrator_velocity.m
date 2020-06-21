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
function [xt,yt] = Integrator_velocity(x0,y0,tspan,options)
Np = numel(x0);               % number of particles
x0 = x0(:); y0 = y0(:);
%% Computing the final positions of the Lagrangian particles:

[~,F] = ode45(@ODEfun,tspan,[x0;y0],options);

xt = F(:,1:end/2);
yt = F(:,end/2+1:end);

    function dy = ODEfun(t,y)
        
        % load the grid over which the velocity is saved
        xi = linspace(0,2*pi,1025);
        yi = linspace(0,2*pi,1025);
        
        N=round(length(y)/2);
        y(1:N,1)     = wrapTo2Pi(y(1:N,1));
        y(N+1:2*N,1) = wrapTo2Pi(y(N+1:2*N,1));
        
        % Interpolate velocity in time and space
        [u1_vec, u2_vec]= interp_vel(t,y,xi,yi);
        dy = zeros(2*N,1);    % a column vector
        dy(1:N,1) = u1_vec(1:N,1);
        dy(N+1:2*N,1) = u2_vec(1:N,1);
        
        function [u_vec, v_vec]=interp_vel(t,y,xi,yi)
            N=round(length(y)/2);
            % load velocity data
            k=floor(t/0.2);
            [ui, vi]=read_vel(k);
            [uf, vf]=read_vel(k+1);
            
            %linear interpolation in time
            u_t = ((k+1)*0.2-t)/0.2*ui + (t-k*0.2)/0.2*uf;
            v_t = ((k+1)*0.2-t)/0.2*vi + (t-k*0.2)/0.2*vf;
            
            %spline interpolation in space
            %u_vec=interp2(x,y, u_t, Y(1:N,1), Y(N+1:2*N,1),’*spline’);
            %v_vec=interp2(x,y, v_t, Y(1:N,1), Y(N+1:2*N,1),’*spline’);
            u_interp = griddedInterpolant({xi,yi},u_t,'spline','none');
            v_interp = griddedInterpolant({xi,yi},v_t,'spline','none');
            u_vec = u_interp(y(1:N,1),y(N+1:2*N,1));
            v_vec = v_interp(y(1:N,1),y(N+1:2*N,1));
            
            function [v1, v2]=read_vel(k)
                
                str1 = '../../Data/turb_u_';
                str2 = pad(int2str(k),4,'left','0');
                str = strcat(str1,str2);
                load(str);
                [n1, n2]=size(u1);
                v1=zeros(n1+1,n2+1); v2=zeros(n1+1,n2+1);
                v1 = [u1 u1(:,1)]; v1=[v1; v1(1,:)]';
                v2 = [u2 u2(:,1)]; v2=[v2; v2(1,:)]';
