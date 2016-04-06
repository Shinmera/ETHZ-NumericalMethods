% Create the matrix C

function C = buildC(A)

n = size(A);
I = eye(n);
C = kron(A,I) + kron(I,A);
