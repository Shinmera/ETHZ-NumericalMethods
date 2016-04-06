function [x_ref,y_ref,T_ref] = refinemesh(x,y,T)
    % MATLAB function for regular refinement of a triangular mesh
    % Arguments:
    % 'x': column vector of x-coordinates of nodes of mesh
    % 'y': column vector of y-coordinates of nodes of mesh
    % 'T': N x 3 - matrix of numbers of nodes of triangles
    % (see documentation of MATLAB's triplot function)
    % Return values:
    % 'x\_ref', 'y\_ref' proivide the coordinates of the nodes
    % of the refined mesh.
    % 'T\_ref' gives the triangles of the refined mesh.

    % First obtain information on the (boundary) edges
    % E and Eb are matrices whose rows contain the numbers
    % of the endpoints of edges
    [E,Eb] = processmesh(T);
    % Append midpoints of edges to list of vertices
    x_ref = [x;0.5*(x(E(:,1))+x(E(:,2)))];
    y_ref = [y;0.5*(y(E(:,1))+y(E(:,2)))];

    % Fetch information about the edges of triangles
    % ET is a matrix whose rows contain the numbers of 
    % the edges of the triangles of the mesh
    ET = getinfo(T,E);
    % Build new list of triangles
    Nt = size(T,1); % Number of triangles
    Nv = length(x); % Number of vertices
    T_ref = zeros(Nt*4,3);
    for l = 1:Nt
      % First son triangle
      T_ref(4*l-3,:) = [T(l,1), ET(l,2)+Nv, ET(l,3)+Nv];
      % Second son triangle
      T_ref(4*l-2,:) = [T(l,2), ET(l,1)+Nv, ET(l,3)+Nv];
      % Third son triangle
      T_ref(4*l-1,:) = [T(l,3), ET(l,1)+Nv, ET(l,2)+Nv];
      % Fourth (interior) son triangle
      T_ref(4*l,:) = [ET(l,1)+Nv, ET(l,2)+Nv, ET(l,3)+Nv];
    end
end
