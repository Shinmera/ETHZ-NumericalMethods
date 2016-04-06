% MATLAB test code for homework problem on "order of convergence from error
% recursion", outputs a table
e(1) = 1; e(2) = 0.8; % Any values <= 1 possible, also random
for k=2:20,  e(k+1) = e(k)*sqrt(e(k-1)); end
le = log(e); diff(le(2:end))./diff(le(1:end-1)),

