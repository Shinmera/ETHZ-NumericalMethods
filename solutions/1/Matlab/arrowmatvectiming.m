nruns = 3; res = [];
ns = 2.^(2:12);
for n = ns
  a = rand(n,1); d = rand(n,1); x = rand(n,1);
  t = realmax;
  for k=1:nruns
    tic; y = arrowmatvec(d,a,x); t = min(toc,t); 
  end;
  res = [res; t];
end
figure('name','timings arrowmatvec');
c1 = sum(res)/sum(ns.^3);
loglog(ns, res,'r+', ns, c1*ns.^3, 'k-');
xlabel('{\bf vector size n}','fontsize',14);
ylabel('{\bf time[s]}','fontsize',14);
title('{\bf timings for arrowmatvec}','fontsize',14);
legend('tic-toc time','O(n^3)','location','best');
print -depsc2 '../PICTURES/arrowmatvectiming.eps';
print -djpeg95 '../PICTURES/arrowmatvectiming.jpg';
