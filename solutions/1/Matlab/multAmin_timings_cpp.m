ps = 4:12;
ns = 2.^ps;
ps2 = 4:13;
ns2 = 2.^ps2;
ts = 1e6*ones(length(ns),3);    % timings
nruns = 10;                     % to average the runtimes

% loop over different Problem Sizes
for j=1:length(ns)
    n = ns(j);
    fprintf('Vector length: %d \n', n);
    x = rand(n,1);
    
    % timing for naive multiplication
    tic;
    for k=1:nruns
      y = min(ones(n,1)*(1:n), (1:n)'*ones(1,n)) * x;
    end
    ts(j,1) = toc/nruns;
    
    % timing multAmin
    tic;
    for k=1:nruns
      y = multAmin(x);
    end
    ts(j,2) = toc/nruns;
    
    % timing multAmin2
    tic;
    for k=1:nruns
      y = multAmin2(x);
    end
    ts(j,3) = toc/nruns;
end

rtime_slow = [3146 3593 7819 32174 176373 560472 1091942 2632146 8337226 42171646] * 1.0e-9;
rtime_slow_loops = [2070 2910 5796 41813 199624 616700 1936121 7330995 26970856 103782319] * 1.0e-9;
rtime_fast = [653 586 639 1070 2156 2893 3648 4879 7365 12883] * 1.0e-9;

c1 = sum(ts(:,2)) / sum(ns);
c2 = sum(ts(:,1)) / sum(ns.^2);

loglog(ns, ts(:,1), '-db', ns, ts(:,2), '-or', ns, ts(:,3), '-xr',...
       ns2, rtime_slow, '-dg', ns2, rtime_slow_loops, '-og', ns2, rtime_fast, '-xy', ...
       ns2, c1*ns2, '-.k', ns2, c2*ns2.^2, '--k', 'linewidth', 2);
legend('Matlab naive','Matlab multAmin','Matlab multAmin2', ...
       'Cpp naive','Cpp naive (with loops)','Cpp fast',...
       'O(n)','O(n^2)', ...
       'Location','NorthWest')

title(sprintf('tic-toc timing averaged over %d runs', nruns),'fontsize',14);
xlabel('{\bf problem size n}','fontsize',14);
ylabel('{\bf runtime (s)}','fontsize',14);

print -depsc2 '../PICTURES/multAmin_timings_cpp.eps';