function p = intpolyval(t,y,x)
% \texttt{t}: row vector of nodes \Blue{$t_{0},\ldots,t_n$}
% \texttt{y}: row vector of data \Blue{$y_{0},\ldots,y_{n}$}
% \texttt{x}: row vector of evaluation points \Blue{$x_{1},\ldots,x_{N}$}
n = length(t); % number of interpolation nodes = degree of polynomial $-1$
N = length(x); % Number of evaluation points stored in \texttt{x}
% Precompute the weights \Blue{$\lambda_{i}$} with effort \Blue{$O(n^{2})$}
for k = 1:n
  lambda(k) = 1 / prod(t(k) - t([1:k-1,k+1:n])); end;
for i = 1:N
% Compute quotient of weighted sums of \Blue{$\frac{\lambda_{i}}{t-t_{i}}$}, effort \Blue{$O(n)$}
  z = (x(i)-t); j = find(z == 0);
  if (~isempty(j)), p(i) = y(j); % avoid division by zero
  else 
    mu = lambda./z; p(i) = sum(mu.*y)/sum(mu);
  end
end
