% compute converg. rate for interpolation by piecewise linear polyn.
% beta-graded mesh in [0,1], singular f(t) = t^alpha, h-convergence
% for different beta
clear all; close all; tic;
alphas  = [0.5, 0.75, 4/3];
betas   = (1:0.5:5);
nn      = 1:50;                              % # used nodes
NumAlph = length(alphas);
NumBeta = length(betas);
NumN    = length(nn);
t_ev    = linspace(0,1,1000);    % points for evaluation and norm
y_ev    = repmat(t_ev(:),1,NumAlph).^repmat(alphas,length(t_ev),1);

Err = zeros(NumAlph, NumBeta, NumN);   % error in max-norm, 3D array
for k=1:NumN;           n = nn(k);
    for m=1:NumBeta;    beta = betas(m);
        t = linspace(0,1,n+1).^beta;              % nodes
        yy = repmat(t(:),1,NumAlph).^repmat(alphas,n+1,1);
        for j=1:NumAlph
            P = PWlineIntp (t,yy(:,j),t_ev);      % interpolation
            [Err(j,m,k), PosErr] = max(abs(y_ev(:,j)-P));
            % PosErr is the index of the point in t_ev with max error
            % LocErr is the index of the subinterval with max error
            LocErr(j,m,k) = sum((t-t_ev(PosErr))<0);
        end
    end
end

rates = zeros(NumAlph,NumBeta);
for j=1:NumAlph
    alpha = alphas(j);
    figure;
    loglog(nn,reshape(Err(j,:,:), NumBeta, NumN),'.-','linewidth',2);
    hold on; %set(gca,'xtick',[3 5 10 20 50])
    title(sprintf('Piecewise linear intp. on graded meshes, alpha =%4.2f',alpha));
    xlabel('n = # subintervals');
    for m=1:NumBeta, leg{m} = sprintf('beta=%4.2f',betas(m));end;
    legend(leg, 'location','nwo');
    
    % estimate the convergence rate, skip the first few n's
    for m =1:NumBeta
        skip = 4;
        ord = polyfit(log(nn(skip+1:end)),...
            reshape(log(Err(j,m,skip+1:end)),1,NumN-skip), 1);
        rates(j,m) = -ord(1);
    end
    
    % plot the convergence rates, the maximum is in re
    axes('Position', [.1, .07, .4, .4]);
    [m1,m2] = max(rates(j,:));
    plot(betas, rates(j,:),'-o',betas(m2),m1,'rs','linewidth',2);
    axis([betas(1),betas(end),0,2.5]);
    xlabel('beta');  ylabel('conv. rate');
    epsname = sprintf('PWlineGraded_%4.2f.eps',alpha);
    print ('-depsc2', epsname); 
end
% display conv. rates
alphas_betas_rates = [0, alphas; betas', rates']
% LocErr
toc