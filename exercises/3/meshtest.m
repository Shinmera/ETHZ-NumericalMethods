%% Test of functions processmesh
figure('name','Triangular mesh');
% Initialize basic mesh data and plot mesh. Note that the variables of
% a MATLAB script are persistent!
meshplotdemo;

% Number of triangles of the mesh
N = numel(x);
M = size(T,1);        

% Obtain relevant information calling a function whose purpose has to
% be determined, see sub-problem 1.
[E,Eb] = processmesh(T);

% Add index numbers to nodes and triangles of the mesh. This is useful
% for understanding the coding of mesh information and also the
% meaning of the data sotred in the matrices E and Eb.
for l=N
  text(x(l),y(l),num2str(l),'fontsize',14,'color','m');
end
for l=1:M
  px = sum(x(T(l,:)))/3;
  py = sum(y(T(l,:)))/3;
  text(px,py,num2str(l),'color','b');
end

% Guess what is plotted and annotated here. This offers a key hint at
% the output of 'processmesh'
k = 1;
for e = E'
  % Note the use of sub-vector extraction in MATLAB.
  plot(x(e),y(e),'k-'); 
  text(0.5*sum(x(e)),0.5*sum(y(e)),num2str(k),'color','k');
  k = k+1;
end
% What does this line of code do? 
for e = Eb', plot(x(e),y(e),'r-'); end

% Obtain some more interesting information on the mesh
disp('Output of getinfo');
ET = getinfo(T,E),

% Loop: Successive smoothing and refinement of the mesh plus graphical
% rendering, which is displayed as Figure 4 of the project sheet.
for l=1:3
  %% refine the mesh and draw it
  [x,y,T] = refinemesh(x,y,T);  
  figure('name',['Refined mesh level ' num2str(l)]);
  triplot(T,x,y,'b-'); title(['Refined mesh level ' num2str(l)]);
  axis([0 1 0 1]); axis equal; axis off;
  hold on; plot(x,y,'r+');
  print('-depsc2',['rmesh' num2str(l) '.eps']);

  %% Smooth the mesh and draw it again
  [x,y] = smoothmesh(x,y,T);
  figure('name',['smoothed mesh level ' num2str(l)]);
  triplot(T,x,y,'b-'); title(['Smoothed mesh level ' num2str(l)]);
  axis([0 1 0 1]); axis equal; axis off;
  hold on; plot(x,y,'r+');
  % print('-depsc2',['smesh' num2str(l) '.eps']);
end
