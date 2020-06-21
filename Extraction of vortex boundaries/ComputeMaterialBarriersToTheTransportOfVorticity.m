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

%% Define parameters of the domain to extract barriers for
xmax = 6.28; ymax = 6.28;
numX = 1100; numY = 1100;

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

%% Load advected positions of initial grid
load('../../particles1100.mat')

% Specify final time to analyze: [t0,t1] = [0,tc]
tc = 49;

%%

load('../../Data/turb_w_0000.mat','w')
w = [w w(:,1)];
w = [w; w(1,:)]';

omega0 = griddedInterpolant({xspan_grid,yspan_grid},w,'spline','none');

str1 = '../../Data/turb_w_';
u_int = 5*tc;
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

%% Extract convex contours larger than a threshold

xspan = xspan_over(1:end-1);
yspan = yspan_over(1:end-1);

Nct = 30;
DeficiencyThresh = 0.05;
MinLength = 0.4;
h1 = contourc(xspan,yspan,vIncrement',Nct);
s = getcontourlines(h1);
bnd = struct('xc',[],'yc',[],'cval',[]);
for k=1:numel(s)
    % check if the contour is closed
    if ( s(k).x(1)==s(k).x(end) ) && ( s(k).y(1)==s(k).y(end) ) 
        Length = sum( sqrt(diff(s(k).x,1).^2+diff(s(k).y,1).^2) );
        if Length>MinLength
            % Check if the contour is convex
             if IsContourConvex(s(k).x,s(k).y,DeficiencyThresh) 
                bnd.xc{end+1} = s(k).x';
                bnd.yc{end+1} = s(k).y'; 
                bnd.cval(end+1,1) = s(k).v';
             end
        end
    end
end

%% Find local maxima of the Lagrangian vorticity-change function
[xp,yp,valp] = Find2DPeak(abs(vIncrement)',xspan,yspan,'maxima');
%- Keeping closed contours that encompass only "1" local maximum.
InMax = cellfun(@(x,y) inpolygon(xp,yp,x,y),bnd.xc,bnd.yc,'UniformoutPut',false);
InMax = cat(1,InMax{:});       % rows --> closed contours & columns --> local maxima
N_InMax = sum(InMax,2);        % Number of local maxima in each contour
indEliminate_1 = N_InMax~=1;
bnd.xc(indEliminate_1)   = [];
bnd.yc(indEliminate_1)   = [];
bnd.cval(indEliminate_1) = [];

[indmax,~] = find( InMax(~indEliminate_1,:)'==1 );
bnd.xp = xp(indmax);
bnd.yp = yp(indmax);
bnd.valp = valp(indmax);

%% Select the outermost contour for each nested family of contours
%- selecting the first point of each contour as a query point
Nct_filt2 = numel(bnd.xc);
if Nct_filt2~=0
    xq = cellfun(@(x) x(1),bnd.xc,'UniformoutPut',true);  
    yq = cellfun(@(x) x(1),bnd.yc,'UniformoutPut',true);  
end

indEliminate_2 = false(Nct_filt2,1);
for ii=1:Nct_filt2
    in = inpolygon(xq,yq,bnd.xc{ii},bnd.yc{ii});  in(ii)=0;
    indEliminate_2(in) = 1;
end
bnd.xc(indEliminate_2)   = [];
bnd.yc(indEliminate_2)   = [];
bnd.cval(indEliminate_2) = [];

%% Plot extracted barriers against the norm of vorticity

str1 = '../../Data/turb_w_';
str2 = pad(int2str(0),4,'left','0');
str = strcat(str1,str2);
load(str);

w = [w w(:,1)]; w = [w; w(1,:)];

h = figure();
AxthicksFnt = 28;
fontsizeaxlab = 28;

pcolor(xspan_grid,yspan_grid,w);
shading interp

hold on

for kk=1:numel(bnd.xc); hold on; plot3(bnd.xc{kk},bnd.yc{kk},1*ones(numel(bnd.xc{kk}),1),'Color',[1 0.1 0.1],'LineWidth',1.8); end

set(gca,'YDir','normal')
set(gcf,'color','w');
set(gcf, 'Position', [0, 0, 845, 845])
axis equal
axis([0 2*pi 0 2*pi])
xlabel('$$x$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
ylabel('$$y$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);

hhF=colorbar(gca);
hhF.Location='eastOutside';
set(get(hhF,'xlabel'),'string','$$\omega$$({\boldmath$${\mathrm{x}}$$}$$,0)$$','Interpreter','latex','FontWeight','normal','FontSize',fontsizeaxlab);
set(hhF,'TickLabelInterpreter', 'latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',AxthicksFnt,'fontWeight','normal');
box on

hold off

