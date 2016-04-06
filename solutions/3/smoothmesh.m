function [xs,ys] = smoothmesh(x,y,T)
    % MATLAB function for smoothing a planar triangular mesh by moving all
    % vertices to the barycenter of their neighboring nodes. The arguments
    % x and y contain the node coordinates, whereas T encodes the triangle
    % node incidence relationship.
    % 
    % Arguments are the same as for 'refinemesh'. The function returns
    % the  new node positions of the smoothed mesh. Note that the
    % connectivity is not affected by the smoothing. 

    % Number of nodes of the mesh
    Nv = length(x); 

    % First obtain the edges of the mesh. E will be a Ne x 2-matrix whose
    % rows give the indices of the nodes adjacent to an edge. Eb contains
    % only the boundary edges in the same format.

    [E,Eb] = processmesh(T); 
    % Number of edges of the mesh
    Ne = size(E,1); 
    % Obtain indices of interior vertices. This is done by first
    % extracting all the vertex numbers that occur in the list of boundary
    % edges. These vertices will be flagged by setting entries in a vector
    % of index numbers to zero.
    idx = (1:Nv)'; idx([Eb(:,1);Eb(:,2)]) = 0; 
    % Indices of vertices in the interior
    int_idx = find(idx > 0); 
    % Indices of vertices on the boundary
    bd_idx = find(idx == 0); 
    % Assemble the (symmetric) graph Laplacian matrix: this is that
    % sparse symmetric quadratic matrix, whose size agrees with the total
    % number of vertices, and whose off-diagonal entries are -1, whenever
    % two vertices are connected by an edge. The diagonal entry is set
    % such that the sum of all entries in each row is equal to zero.
    L = sparse(E(:,1),E(:,2),ones(Ne,1),Nv,Nv); L = L+L';
    L = spdiags(sum(L)',[0],Nv,Nv) - L;
    % Remove rows and columns corresponding to vertices on the boundary
    % using MATLAB's sub-matrix extraction capability.
    A = L(int_idx,int_idx);
    % Build right hand side vector. This vector depends on the positions
    % of the nodes on the boundary. It can be obtained by multiplying the
    % vector of either x or y coordinates of nodes on the boundary with
    % that columns of the graph Laplacian matrix that corresponds to nodes
    % on the boundary.
    B = [x(bd_idx),y(bd_idx)]; 
    B = -L(int_idx,bd_idx)*B;
    % Solve linear system for both x- and y-coordinates simultaneously
    X = A\B;
    % Build output vectors by merging new positions of interior nodes and
    % fixed positions of nodes on the boundary.
    xs = x; ys = y; 
    xs(int_idx) = X(:,1); ys(int_idx) = X(:,2);
end
