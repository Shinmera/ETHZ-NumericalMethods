% compute converg. rate for interpolation by piecewise linear polyn.
% uniform mesh in [0,1], singular f(t) = t^alpha, h-convergence

clear all; close all; tic;
alphas = [0.1:0.2:2.9];
nn = 1:50;                              % # used nodes
NumAlph = length(alphas);
NumN = length(nn);
t_ev = linspace(0,1,1000);    % points for evaluation and norm
y_ev = repmat(t_ev(:),1,NumAlph).^repmat(alphas,length(t_ev),1);

Err = zeros(NumAlph, NumN);             % error in max-norm
LocErr = zeros(NumAlph, NumN);          % location of maximal error
for k=1:NumN
    n = nn(k);
    t = linspace(0,1,n+1);              % nodes
    yy = repmat(t(:),1,NumAlph).^repmat(alphas,n+1,1);
    for j=1:NumAlph
        P = PWlineIntp (t,yy(:,j),t_ev);      % interpolation
        [Err(j,k), PosErr] = max(abs(y_ev(:,j)-P));
        % PosErr is the index of the point in t_ev with max error
        % LocErr is the index of the subinterval with max error
        LocErr(j,k) = sum((t-t_ev(PosErr))<0);  
        % warning if the maximal error is not where expected
        if (alphas(j)<2 && LocErr(j,k)~=1) || (alphas(j)>2 && LocErr(j,k)~=n)
        warning('(alpha=%4.2f, N=%d), max. err. in interval %d', alphas(j),n,LocErr(j,k));
        end 
    end
end

figure;  
loglog(nn,Err,'.-','linewidth',2);
hold on; set(gca,'xtick',[3 5 10 20 50])
title('Piecewise linear intp. on uniform meshes: error in max-norm');
xlabel('n = # subintervals');
for j=1:NumAlph, leg{j} = sprintf('alpha=%4.2f',alphas(j)); end;
legend(leg, 'location','nwo');

% estimate the convergence rate
rates = zeros(1,NumAlph);
for j =1:NumAlph
    ord = polyfit(log(nn), log(Err(j,:)), 1);
    % check the use of polyfit:
    %loglog([1,5], [1,5^ord(1)]*exp(ord(2))/2,'k--' );
    rates(j) = -ord(1);
end

% plot the convergence rate
axes('Position', [.1, .1, .3, .2]);
plot(alphas, rates,'-o','linewidth',2);
axis([alphas(1),alphas(end),0,2.2]);
xlabel('alpha');  ylabel('conv. rate');
hold on; plot([0 2],[0 2],'r', [2,alphas(end)],[2,2],'r');
print -depsc2 'PWlineConv.eps';

alpha_orders = [alphas(:),rates(:)]  % display conv. rates
toc