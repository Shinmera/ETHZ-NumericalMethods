function y = Kron_B(A,B,x)
% Return y = kron(A,B)*x  (smart version)
% Input:   A,B:     2 n x n matrices.
%          x:       Vector of length n*n.
% Output:  y:       Result vector of length n*n.

% check size of A
[n,m] = size(A);
assert(n == m, 'expected quadratic matrix')

% kron gives a matrix with n x n blocks
% block i,j is A(i,j)*B
% => y = M*x can be done block-wise so that we reuse B*x(...)

% init
y = zeros(n*n,1);
% loop first over columns and then (!) over rows
for j = 1:n
    % reuse B*x(...) part (constant in given column) => O(n^2)
    z = B*x((j-1)*n+1:j*n);
    % add to result vector (need to go through full vector) => O(n^2)
    for i = 1:n
        y((i-1)*n+1:i*n) = y((i-1)*n+1:i*n) + A(i,j)*z;
    end
end
% Note: complexity is O(n^3)
end
