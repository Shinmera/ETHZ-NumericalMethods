function X = Pow (A, k)
% Pow - Return A^k for square matrix A and integer k
% Input:    A:    n*n matrix
%           k:    positive integer
% Output:   X:    n*n matrix X = A^k

% transform k in basis 2 
bin_k = de2bi(k) ;
M = length(bin_k);
X = eye(size(A));

for j = 1:M 
    if bin_k(j) == 1
        X = X*A;
    end
    A = A*A;    % now A{new} = A{initial} ^(2^j)
end
