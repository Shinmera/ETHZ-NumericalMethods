function GaussConv(f_hd)
if nargin<1;    f_hd = @(t) sinh(t);    end;
I_exact = quad(@(t) asin(t).*f_hd(t),-1,1,eps)

n_max = 50;  nn = 1:n_max;
err = zeros(size(nn));
for j = 1:n_max
    [x,w] = gaussquad(nn(j));
    I = dot(w, asin(x).*f_hd(x));
    % I = GaussArcSin(f_hd,nn(j));   % using pcode
    err(j) = abs(I - I_exact);
end

close all; figure; 
loglog(nn,err,[1,n_max],[1,n_max].^(-3),'--','linewidth',2);    
title('{\bf Convergence of Gauss quadrature}','fontsize',14);
xlabel('{\bf n = # of evaluation points}','fontsize',14);
ylabel('{\bf error}','fontsize',14);
legend('Gauss quad.','O(n^{-3})');
print -depsc2 '../PICTURES/GaussConv.eps';
