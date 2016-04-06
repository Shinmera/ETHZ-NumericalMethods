clear all; close all;
nruns = 10;                % we average over a few runs
nn = 30:3:60;              % dimensions used
kk = (2:50);               % powers used
tt1 = zeros(length(nn), length(kk));  % times for Pow
tt2 = zeros(length(nn), length(kk));  % times for Matlab power

for i = 1:length(nn)
    n = nn(i);
    
    % matrix with |eigenvalues| =1:
    % A = vander([exp(1i * 2 * pi * [1:n]/n)])/sqrt(n);
    A = exp(1i * 2 * pi * [1:n]'*[1:n]/n)/sqrt(n);
    
    for j=1:length(kk)
        k = kk(j);
        tic
        for run = 1:nruns
            X = Pow(A, k);
        end
        tt1(i,j) = toc;
        
        tic
        for run = 1:nruns
            XX = A^k;
        end
        tt2(i,j) = toc;
        n_k_err = [n, k, max(max(abs(X-XX)))]
    end
    
end

figure('name','Pow timings');
subplot(2,1,1)
n_sel=6;             %plot in k only for a selected n
% expected logarithmic dep. on k, semilogX used:
semilogx(kk,tt1(n_sel,:),'m+',     kk,tt2(n_sel,:),'ro',...
       kk,sum(tt1(n_sel,:))*log(kk)/(length(kk)*log(k)), 'linewidth', 2);
xlabel('{\bf power k}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
title(sprintf('tic-toc timing averaged over %d runs, matrix size = %d',...
     nruns, nn(n_sel)),'fontsize',14);
legend('our implementation','Matlab built-in',...
       'O(C log(k))','location','northwest');

subplot(2,1,2)
k_sel = 35;      %plot in n only for a selected k
loglog(nn, tt1(:,k_sel),'m+',     nn, tt2(:,k_sel),'ro',...
       nn, sum(tt1(:,k_sel))* nn.^3/sum(kk.^3), 'linewidth', 2);
xlabel('{\bf dimension n}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
title(sprintf('tic-toc timing averaged over %d runs, power = %d',...
     nruns, kk(k_sel)),'fontsize',14);
legend('our implementation','Matlab built-in',...
       'O(n^3)','location','northwest');
print -depsc2 'Pow_timings.eps';
