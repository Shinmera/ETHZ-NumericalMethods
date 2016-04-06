function y = multAmin2(x)
% O(n), no-for version
n = length(x);
v = cumsum(x(end:-1:1));
w = cumsum(x.*(1:n)');
y = w + (1:n)'.*[v(n-1:-1:1);0];
