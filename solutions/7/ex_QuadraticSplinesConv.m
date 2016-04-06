function ex_QuadraticSplinesConv(f,n_mesh,N)
tic
if nargin<3
    f = @(t) exp(sin(2*pi*t));   % function to be interpolated
    n_mesh = 10;                 % number of meshes used
    N = 10000;                   % number of evaluation points
end
nn = 2.^(2:1+n_mesh);
x = linspace(0,1,N);     % N evaluation points
x = x(0 <= x & x < 1);   % be sure x belongs to [0,1):
f_x = f(x);              % exact values of f to compute error
nruns = 3;
    
err = zeros(1, n_mesh);
timings = realmax*ones(1, n_mesh);

for j=1:n_mesh
    n = nn(j); 
    t = 1/n:1/n:1-1/n;   % n-1 equispaced nodes in (0,1)
    middle_points = ([t,1]+[0,t])/2;
    y = f(middle_points);
    
    for r = 1:nruns
        tic;
        % compute and evaluate the quadratic spline
        eval_x = quadspline(t,y,x);
        timings(j) = min(toc, timings(j));
    end
    err(j) = max(abs(eval_x - f_x));
end
close all; figure; subplot(1,2,1);
loglog(nn, err,'o-', nn,nn.^(-3), 'k',  'linewidth',2);
legend('Error', 'O( n^{-3} )');grid on;
xlabel('Number of intervals','fontsize',14);
ylabel('Error in maximum norm','fontsize',14);
subplot(1,2,2);
loglog(nn, timings,'o-', nn,nn/10000, 'k','linewidth',2);
legend('Time','O(n)'); grid on;
xlabel('Number of intervals','fontsize',14); 
ylabel('Time for computation and evaluation','fontsize',14);
print -depsc2 'time-error_QuadSplines.eps';
max_n = nn(n_mesh)
total_time = toc