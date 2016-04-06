n = 10;              b = 2;
t_eq = (0:n)/n;      t_gr =t_eq.^b;
close all; plot([0:0.01:1],(0:0.01:1).^b,'r','linewidth',2); hold on; 
for j=1:n+1;   plot([t_eq(j),t_eq(j),0], [0,t_gr(j),t_gr(j)],'k');
    plot([t_eq(j),0], [0,t_gr(j)],'ro','linewidth',2);  end;
axis square;  xlabel('uniform mesh','fontsize',14);
ylabel('algeb. graded mesh, beta=2','fontsize',14);
set(gca,'xtick',t_eq,'ytick',t_gr);print -depsc2 'GradedMesh.eps';