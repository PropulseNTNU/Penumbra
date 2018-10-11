
function particle_dynamics_images(x_in,q_in,t_in,tskip,nvect,filename)

global lam d

% Produce multiple images called 'filename-i.png'

%% Plotting options and parameters
% Plot streamlines and or vector field arrows?
plotstream = false ;
plotvect   = true  ;

plotpath   = true  ; % False will plot over the surface

overide_range = true;
xmin = -0.2; 
xmax =  1.2;
ymin = -1.2;
ymax =  0.2;
zmin =  -1.4;
zmax =  0;

% range of angles to sweep through
th1i = 32;
th1f = 64;
th2i = 0;
th2f = 20;

n = 5 ; % number of streamlines/velocity arrows to plot over whole domain (not for path vectors)

psize = 0.06; % size of particle relative to domain

if overide_range
    range = [xmin xmax ymin ymax zmin zmax];
end

%% Position euler angles from input
nT = floor(length(t_in)/tskip); % scaled length of solutions
xscale = x_in(1:tskip:end-tskip,:);
qscale = q_in(1:tskip:end-tskip,:);
tscale = t_in(1:tskip:end-tskip);

x = xscale(:,1); y = xscale(:,2); z = xscale(:,3);

ea = zeros(nT,3);

if tskip==1 % Fixes bug that doesnt allow for tskip=1
    nT=nT-1;
end

for it=1:nT
    ea(it,:) = 180/pi*quat2ea(qscale(it,:)');
end

%% Axis
xrange = [min(x) max(x)]; yrange = [min(y) max(y)]; zrange = [min(z) max(z)];

% Length of each axis
xlen = max(xrange)-min(xrange); ylen = max(yrange)-min(yrange); zlen = max(zrange)-min(zrange);

% find the longest axis then make all other axis equal (cubic domain)
if xlen>ylen && xlen>zlen %if x axis is the longest, then scale other two to have equal length
    xy=(xlen-ylen)/2;
    yrange = [min(y)-xy max(y)+xy];
    xz=(xlen-zlen)/2;
    zrange = [min(z)-xz max(z)+xz];
elseif ylen>xlen && ylen>zlen % if y is longest then do same...
    yx=(ylen-xlen)/2;
    xrange = [min(x)-yx max(x)+yx];
    yz=(ylen-zlen)/2;
    zrange = [min(z)-yz max(z)+yz];
elseif zlen>xlen && zlen>ylen
    zx=(zlen-xlen)/2;
    xrange = [min(x)-zx max(x)+zx];
    zy=(zlen-ylen)/2;
    yrange = [min(y)-zy max(y)+zy];
end
ax = max(xrange)-min(xrange); % Length of each axis

if overide_range
    ax = range(2)-range(1);
end

if ~ overide_range
    range = [xrange yrange zrange]; % Plotting range
end
%% Particle dimensions

% Relative size of particle
rpsize = psize*ax;

% Aspect ratios of particle
if lam>1
    arx = 1.5/lam; ary = 1.5/lam; arz = 1.5;
elseif lam<1
    arx = 1; ary = 1; arz = lam;
elseif lam==1
    arx = 1; ary = 1; arz = 1;
end
arx = arx*rpsize; ary = ary*rpsize; arz = arz*rpsize;

%% Velocity field

if plotpath
    % up is the velocity vector along the particle path (with time dependence)
    up = zeros(3,nT);
    for it=1:nT
        if velfield == 0
            [up(:,it),~] = getTurbFluid([x(it),y(it),z(it)],tscale(it));
        else
            [up(:,it),~] = analytic_fluid_field(x(it),y(it),z(it),tscale(it));
        end
    end
    up=up';
    
else
    % vel field for the whole domain at t=0:
    [xx,yy,zz,u,v,w,sx,sy,sz] = plot_velfield(n,range(1),range(2),range(3),range(4),range(5),range(6));
    
end
%% Animate particle
for it=1:nT
    
    % Define ellipsoid at origin
    [xs, ys, zs] = ellipsoid(0,0,0,arx,ary,arz,10);
    S = surfl(xs, ys, zs);
    s.EdgeColor = 'none';
    s.FaceColor = 'flat';
    
    % Rotate particle using Euler angles
    rotate(S,[1 0 0],real(ea(it,1)));
    rotate(S,[0 1 0],real(ea(it,2)));
    rotate(S,[0 0 1],real(ea(it,3)));
    
    xs = S.XData; ys = S.YData; zs = S.ZData;
    
    % Plot and translate particle to correct position
    surfl(xs+x(it), ys+y(it), zs+z(it));
    title(['Frame ',num2str(it),' out of ',num2str(nT)])
    xlabel('x'); ylabel('y'); zlabel('z');
    colormap jet
    
    % Plot velocity feild/streamlines and path
    hold on
    if plotvect
        if plotpath % plot vecotr field along path:
            quiver3(x(1:nvect:end),y(1:nvect:end),z(1:nvect:end),up(1:nvect:end,1),up(1:nvect:end,2),up(1:nvect:end,3),'black','AutoScale','on','AutoScaleFactor',2*exp(-d^2*tscale(it)))
        else  % plot vector field on surface of domain:
            quiver3(xx,yy,zz,u,v,w,'black','AutoScale','on', 'AutoScaleFactor',exp(-d^2*tscale(it)))
        end
    end
        
    % Plot the path
    plot3(x,y,z,'k-','linewidth',2)
    view([th1i+it*(th1f-th1i)/nT,th2i+it*(th2f-th2i)/nT])
    
    hold off
    axis(range)
    set(gca,'fontsize',16)

    pause(0.001)
    
    filename_i = sprintf('%s-%i',filename,it);
    print(filename_i,'-dpng')
    
end