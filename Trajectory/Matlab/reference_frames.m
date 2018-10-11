% reference_frames
close all
addpath(genpath('../functions'));

% particle location:
xp = 0.; yp = 0.3; zp = 0.3;

% particle euler angles
ea = [-pi/8,-pi/8,0];

Q = ea2mat(ea);

relsize = 2; % relative size of stationary frame axis

% Translating frame:
T = [ 1  0  0;
      0  1  0;
      0  0  1];
  
% Stationary frame:
S = 1*[ 1  0  0;
      0  1  0;
      0  0  1];

% Body frame:
B=T;
for i=1:3
    B(i,:)=(Q*T(i,:)')';
end

x = xp*ones(3,1);
y = yp*ones(3,1);
z = zp*ones(3,1);

o = zeros(3,1);
psize = 1.3;

[xs, ys, zs] = ellipsoid(xp,yp,zp,psize*0.03,psize*0.03,psize*0.1,20);
P = surfl(xs, ys, zs);
rotate(P,[1 0 0],real(180/pi*ea(1)));
rotate(P,[0 1 0],real(180/pi*ea(2)));
rotate(P,[0 0 1],real(180/pi*ea(3)));
xs = P.XData; ys = P.YData; zs = P.ZData;

close all


f=figure(1);
hold on
quiver3(o,o,o,S(:,1),S(:,2),S(:,3),'AutoScaleFactor',1,'color','k','linewidth',4)
quiver3(x,y,z,T(:,1),T(:,2),T(:,3),'AutoScaleFactor',1 ,'color','k','linewidth',2,'lineStyle','-')
quiver3(x,y,z,B(:,1),B(:,2),B(:,3),'AutoScaleFactor',1 ,'color','b','linewidth',2,'lineStyle','-')
P = surf(xs, ys, zs);
P.EdgeColor = 'none';
P.FaceColor = 'flat';
% legend('\bf{x}_s','\bf{x}_t','\bf{x}')
view([115 30])
hold off
grid on 
set(gca,'XTick',[]); % hides the ticks on x-axis
set(gca,'color','None'); % hides the white bckgrnd
set(gca,'YColor','w'); % changes the color of y-axis to white
set(gca,'visible','off'); %hides all axes, ticks and labels
set(findall(gca, 'type', 'text'), 'visible', 'on'); % Keeps the title if we need one

xs = 0.1; ys = 0; zs = 1; txts = '\bf{x}_s';
txs = text(xs,ys,zs,txts,'HorizontalAlignment','right');

xt = xp; yt = yp; zt = zp; txtt = '\bf{x}_t';
txt = text(xt+0.1,yt,zt+1,txtt,'HorizontalAlignment','right');

xxb = Q*[0.15, 0, 1]';
xb = xxb(1); yb = xxb(2); zb = xxb(3);

txtb = '\bf{x}';
txb = text(xb+xp,yb+yp,zb+zp,txtb,'HorizontalAlignment','right');

txs.FontSize = 20;
txt.FontSize = 20;
txb.FontSize = 20;


% printpdf(f,'/Users/benjakt/Box Sync/PhD/Figures/reference_frames')



% saveas(f,'/Users/benjakt/Box Sync/PhD/Figures/reference_frames','eps')
