nruns = 3; res = [];
ns = 2.^(2:12);
for n = ns
  a  = rand(n,1); d = rand(n,1); x = rand(n,1);
  t  = realmax;
  t2 = realmax;
  for k=1:nruns
    tic; y  = arrowmatvec(d,a,x);  t  = min(toc,t);
    tic; y2 = arrowmatvec2(d,a,x); t2 = min(toc,t2);
  end;
  res = [res; t t2];
end
figure('name','timings arrowmatvec and arrowmatvec2');
c1 = sum(res(:,1))/sum(ns.^3);
c2 = sum(res(:,2))/sum(ns);
loglog(ns, res(:,1),'r+', ns, res(:,2),'bo',...
       ns, c1*ns.^3, 'k-', ns, c2*ns, 'g-');
xlabel('{\bf vector size n}','fontsize',14);
ylabel('{\bf time[s]}','fontsize',14);
title('{\bf timings for arrowmatvec and arrowmatvec2}','fontsize',14);
legend('arrowmatvec','arrowmatvec2','O(n^3)','O(n)',...
       'location','best');
print -depsc2 '../PICTURES/arrowmatvec2timing.eps';
