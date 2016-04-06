function ex_QuadraticSplinesPlot(f,n,N, func)

if nargin<3
    f = @(t) exp(sin(2*pi*t));   % function to be interpolated
    n = 10;                      % number of subintervals
    N = 200;                     % number of evaluation points
    func = @quadspline_better;
end

t = 1/n:1/n:1-1/n;               % n-1 equispaced nodes in (0,1)
x = linspace(0,1,N);             % N evaluation points
x = x(0 <= x & x < 1);           % be sure x belongs to [0,1)
middle_points = ([t,1]+[0,t])/2; % computes the n middle points [0,t,1]
y = f(middle_points);            % evaluate f on the middle_points

% compute and evaluate the quadratic spline
eval_x = func(t,y,x);
%eval_x_better = quadspline_better(t,y,x);

close all; figure
plot(middle_points,y,'ro',x,eval_x,'k',x,f(x),'--','linewidth',2);
legend('Data','Quadratic spline','Function f');
print -depsc2 'interpol_quad_T2.eps'

end

% test:
% ex_QuadraticSplinesPlot(@(t)sin(2*pi*t).^3, 10,500)
