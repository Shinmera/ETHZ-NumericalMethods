ps = 4:12;
ns = 2.^ps;
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

c1 = sum(ts(:,2)) / sum(ns);
c2 = sum(ts(:,1)) / sum(ns.^2);

loglog(ns, ts(:,1), '-k', ns, ts(:,2), '-og', ns, ts(:,3), '-xr',...
       ns, c1*ns, '-.b', ns, c2*ns.^2, '--k', 'linewidth', 2);
legend('naive','multAmin','multAmin2','O(n)','O(n^2)',...
       'Location','NorthWest')
title(sprintf('tic-toc timing averaged over %d runs', nruns),'fontsize',14);
xlabel('{\bf problem size n}','fontsize',14);
ylabel('{\bf runtime (s)}','fontsize',14);

print -depsc2 '../PICTURES/multAmin_timings.eps';
