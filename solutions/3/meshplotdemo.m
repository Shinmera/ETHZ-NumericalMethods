% MATLAB demonstration for visualizing a planes triangular mesh
% Initialize node coordinates
% First the x-coordinates
x = [1.0;0.60;0.12;0.81;0.63;0.09;0.27;0.54;0.95;0.96];
% Next the y-coordinates
y = [0.15;0.97;0.95;0.48;0.80;0.14;0.42;0.91;0.79;0.95];
% Then specify triangles through the indices of their vertices. These
% indices refer to the ordering of the coordinates as given in the
% vectors x and y.
T = [8 2 3;6 7 3;5 2 8;7 8 3;7 5 8;7 6 1;...
     4 7 1;9 5 4;4 5 7;9 2 5;10 2 9];
% Call the plotting routine; draw mesh with blue edges
triplot(T,x,y,'b-'); title('A simple planar triangular mesh');
xlabel('{\bf x}'); ylabel('{\bf y}'); 
axis([-0.05 1.05 -0.05 1.05]); axis equal; 
% Mark nodes with red stars
hold on; plot(x,y,'r*');

% Save plot a vector graphics
print -depsc2 'meshplot.eps';
