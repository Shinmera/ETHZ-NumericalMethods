function y = multAmin(x)
% O(n), slow version
n = length(x);
y = zeros(n,1);
v = zeros(n,1);
w = zeros(n,1);

v(1) = x(n);
w(1) = x(1);
for j = 2:n
    v(j) = v(j-1)+x(n+1-j);
    w(j) = w(j-1)+j*x(j);
end
for j = 1:n-1
    y(j) = w(j) + v(n-j)*j;
end
y(n) = w(n);

% To check the code, run:
% n=500; x=randn(n,1); y = multAmin(x);
% norm(y - min(ones(n,1)*(1:n), (1:n)'*ones(1,n)) * x)
