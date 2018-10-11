function [x,y,z,u,v,w,sx,sy,sz] = plot_velfield(m,xmin,xmax,ymin,ymax,zmin,zmax)

global velfield

% generates the coordinates (x,y,z) and the vector field (u,v,w) for 
% the solution to the Navier-Stokes equation. 
% m is the number of streamlines in each dimension, ax is half the length
% of the cubic domain, a and d are paramters that control the field

t = 0;

lx = linspace(xmin,xmax,m);
ly = linspace(ymin,ymax,m);
lz = linspace(zmin,zmax,m);
[x,y,z] = meshgrid(lx,ly,lz);

% u = -a*(exp(a*x).*sin(a*y+d*z)+exp(a*z).*cos(a*x+d*y)).*exp(-d^2*t);
% v = -a*(exp(a*y).*sin(a*z+d*x)+exp(a*x).*cos(a*y+d*z)).*exp(-d^2*t);
% w = -a*(exp(a*z).*sin(a*x+d*y)+exp(a*y).*cos(a*z+d*x)).*exp(-d^2*t);

u = zeros(m); v = zeros(m); w = zeros(m); 
for ix=1:m
    for iy=1:m
        for iz=1:m
            if velfield == 0 
                [U,~] = getTurbFluid([lx(ix),ly(iy),lz(iz)],t);
            else
                [U,~] = analytic_fluid_field(lx(ix),ly(iy),lz(iz),t);
            end
            u(ix,iy,iz) = U(1); v(ix,iy,iz) = U(2); w(ix,iy,iz) = U(3);
        end
    end
end

sx = x; sy = y; sz = z;

end

