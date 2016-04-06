f = @(t) exp(sin(2*pi*t)); 
sizes = 2.^(2:15);
N = 1000000;

times = [];
times_better = [];

for n = sizes,
    t = 1/n:1/n:1-1/n;               % n-1 equispaced nodes in (0,1)
    x = linspace(0,1,N);             % N evaluation points
    x = x(0 <= x & x < 1);           % be sure x belongs to [0,1)
    middle_points = ([t,1]+[0,t])/2; % computes the n middle points [0,t,1]
    y = f(middle_points);            % evaluate f on the middle_points
    
    tic
    quadspline(t,y,x);
    tm = toc;
    times = [times, tm];
    tic
    quadspline_better(t,y,x);
    tm = toc;
    times_better = [times_better, tm];
end

figure
loglog(sizes, times, sizes, times_better,  'linewidth',2);
title('Runtime comparison beween SMW/no-SMW')
xlabel('interpolation nodes')
ylabel('runtime')
legend('No SMW fomrula', 'SMW formula');