function y = Kron_C (A,B,x)
% Return y = kron(A,B)*x  (smart version with reshapes)
% Input:   A,B:       2 n x n matrices.
%            x:       Vector of length n*n.
% Output:    y:       Result vector of length n*n.

% check size of A
[n,m] = size(A);  assert(n == m, 'expected quadratic matrix')

% init
yy = zeros(n,n);
xx = reshape(x,n,n);

% precompute the multiplications of B with the parts of vector x
Z = B*xx;
for j=1:n 
    yy = yy + Z(:,j) * A(:,j)';
end
y = reshape(yy,n^2,1);
end
% Notes: complexity is O(n^3)