function ET = getinfo(T,E)
    % T: N x 3 matrix, nodes of triangles
    % (N = number of triangles of the mesh)
    % E: L x 2 matrix, endpoints of edges
    % (L = number of edges of the mesh)
    % ET: is a N x 3 -matrix (N = no. of triangles)
    % whose rows contain the index numbers of the edges of the 
    % the triangles. 
    L = size(E,1);   % Number of edges
    N = max(max(T)); % Number of nodes
    % Create a sparse N x N - matrix A, for which A(i,j) is zero, if the
    % nodes i and j are not connected by an edge, and for which A(i,j) gives
    % the INDEX of the edge from i to j, if it exists. The index of an edge
    % corresponds to its row number in the E matrix
    A = sparse(E(:,1),E(:,2),(1:L)',N,N); 
    % Loop over all triangles
    ET = zeros(size(T,1),3);
    i = 1;
    for tri=T'
      % The node numbers of the current triangles are contained
      % in the loop variable tri (row vector of length 3)
      % Eloc is a 3x3-matrix that contains the numbers of 
      % the edges of the triangles: Eloc(i,j) gives the number
      % of the edge connecting vertex i and vertex j, 1<=i,j<=3,
      % of the triangle
      Eloc = full(A(tri,tri)); Eloc = Eloc + Eloc';
      % The three edge numbers are stored in a row of the matrix
      % ET: order: edge [2,3] - edge [1,3] - edge [1,2]
      ET(i,:) = Eloc([8 7 4]);
      i = i + 1;
    end
end
