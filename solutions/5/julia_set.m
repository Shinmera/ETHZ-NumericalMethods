function julia_set(N_it, num_grid_points)
close all;
if nargin<2;   N_it = 25;  num_grid_points = 200; end;

% initialize matrix for colors
col = ones(num_grid_points);
a = 2.0;
xx = linspace(-a,+a,num_grid_points);
yy = linspace(-a,+a,num_grid_points);

% roots of unity:
z1 = [ 1;    0];
z2 = [-1/2;  sqrt(3)/2];
z3 = [-1/2; -sqrt(3)/2];

for ny = 1:length(yy)
    for nx = 1:length(xx)
        % for each starting point in the grid
        v = [xx(nx); yy(ny)];
        
        for k=1:N_it;
            F  = [v(1)^3-3*v(1)*v(2)^2-1; 3*v(1)^2*v(2)-v(2)^3];
            DF = [3*v(1)^2-3*v(2)^2, -6*v(1)*v(2);
                  6*v(1)*v(2)      , 3*v(1)^2-3*v(2)^2];
            v = v - DF\F;           % Newton update
            
            if     ((norm(v-z1)) < 1e-4)
                col(ny,nx) = 1 + k;                 break;
            elseif ((norm(v-z2)) < 1e-4)
                col(ny,nx) = 1 + k + N_it;          break;
            elseif ((norm(v-z3)) < 1e-4)
                col(ny,nx) = 1 + k + 2*N_it;        break;
            end
        end
    end
end
st = 1/N_it;         % build a RGB colormap,
% 1s at the beginning to recognize slowly-converging points (white)
mycolormap = [ [1,1-st:-st:0, zeros(1,2*N_it)];
               [1,zeros(1,N_it), 1-st:-st:0, zeros(1,N_it)];
               [1,zeros(1,2*N_it), 1-st:-st:0] ]';
% mycolormap = 1 - mycolormap;    % cmy, pyjamas version...
% mycolormap = [ [1,1-st:-st:0, st:st:1,    st:st:1,];
%                [1,st:st:1,    1-st:-st:0, st:st:1,];
%                [1,st:st:1,    st:st:1,    1-st:-st:0] ]';
colormap(mycolormap);            % built in: colormap jet(256);
% this is the command that creates the plot:
pcolor(col);                     
caxis([1,3*N_it+1]); shading flat; axis('square','equal','off');
print -depsc2 'ex_JuliaSet.eps';
