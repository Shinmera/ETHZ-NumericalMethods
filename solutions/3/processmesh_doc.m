function [E,Eb] = processmesh(T)
    % PROCESSMESH function extracting edge information of a mesh
    % Argument 'T' is a N x 3 matrix containing the vertex numbers of
    % each of the N triangles of a triangular mesh.
    % Return value E is a matrix whose rows correspond to the edges of the
    % mesh and give the indices of the endpoints of the edges.  Return
    % value Eb is the same data structure for edges on the boundary.

    % Number of nodes of the triangular mesh
    N = max(max(T)); 
    % Number of triangles of the mesh 
    M = size(T,1);   
    % Sorts the vertices of each triangle in ascending order, sort by row
    T = sort(T,2);  
    % Create a matrix C with two columns; the first column contains the
    % numbers of first points of edges, the second the number of second
    % points. Each row is an edge
    C = [T(:,1) T(:,2); T(:,2) T(:,3); T(:,1) T(:,3)];
    % The matrix C is used to initialize a sparse matrix A through
    % MATLAB's 'sparse' command. The rows and the columns of the matrix
    % correspond to nodes of the mesh. If A(i,j) is different from zero
    % the two nodes with numbers i and j are connected by an edge. Note
    % that A(i,j) is equal to 2, if [i,j] is an interior edge. For a
    % boundary edge [i,j] A(i,j) == 1, because an interior edge is
    % contained twice in the matrix C.
    A = sparse(C(:,1), C(:,2), ones(3*M,1), N, N);
    % MATLAB's find command gives all pairs of row and column indices of
    % matrix elements satisfying a certain condition. This is used to
    % build the matrices E and Eb
    [I,J] = find(A > 0); E = [I,J];
    [I,J] = find(A == 1); Eb = [I,J];
end
