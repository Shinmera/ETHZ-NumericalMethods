%% Kron (runtimes)
clear all; close all;
nruns = 3;  % we average over a few runs
N = 2.^(2:10);
T = zeros(length(N),5);
for i = 1:length(N)
    n = N(i)    % problem size n
    % % initial matrices A,B and vector x
    % A = reshape(1:n*n, n, n)';
    % B = reshape(n*n+1:2*n*n, n, n)';
    % x = (1:n*n)';    
    % alternative choice:
    A=rand(n);    B=rand(n);    x=rand(n^2,1);
     
    tic;    % smart implementation #1
    for ir = 1:nruns
        yB = Kron_B(A,B,x);
    end
    tb = toc/nruns;
  
    tic;   % smart implementation #2: with reshapes
    for ir = 1:nruns
        yC = Kron_C(A,B,x);
    end
    tc = toc/nruns;
    fprintf('Error B-vs-C:   %g\n', norm(yB-yC))   
    
    tic;   % implementation with kron and matrix*vector
    if(N(i)<128)
        for ir = 1:nruns
            yA = kron(A,B)*x;
        end
        ta = toc/nruns;
        fprintf('Error A-vs-B:   %g\n', norm(yA-yB))
    else
        ta=0;
    end        
        
    tic;  % implementation with 1-line command!
          % inspired by the solution of Kevin Bocksrocker
    for ir = 1:nruns   
        yD = reshape(B*reshape(x,n,n)*A',n^2,1); 
    end
    td = toc/nruns;
    fprintf('Error D-vs-B:   %g\n', norm(yD-yB));   
       
    fprintf('Timings: Matrix: %g\n         Smart:  %g\n',ta,tb)
    fprintf('         Smart+Reshape:  %g\n         1-line Smart:  %g\n',tc,td)
    T(i,:) = [n, ta, tb, tc, td];
end

a = 2.^(2:7);
b = 2.^(2:10);
kron_slow_cpp = [3427 12702 110188 1457921 11898663 179823481] * 1e-9;
kron_fast_cpp = [1831 2999 8331 41265 196635 886109 2791087 19921186 236021413 ] * 1e-9;
kron_super_fast_cpp = [1875 2120 3974 18675 91372 396124 1248023 7067393 45488499] * 1e-9;


% log-scale plot for investigation of asymptotic complexity
a2 = sum(T(:,3)) / sum(T(:,1).^2);
a3 = sum(T(:,3)) / sum(T(:,1).^3);
a4 = sum(T(1:5,2)) / sum(T(1:5,1).^4);
figure('name','kron timings');
loglog(T(:,1),T(:,2),'m+',  T(:,1),T(:,3),'ro',...
       T(:,1),T(:,4),'bd',  T(:,1),T(:,5),'gp',...
       a, kron_slow_cpp, 'm', b, kron_fast_cpp, 'r', b, kron_super_fast_cpp, 'g', ...
       T(:,1),(T(:,1).^2)*a2,'k-',  T(:,1),(T(:,1).^3)*a3,'k--',...
       T(1:5,1),(T(1:5,1).^4)*a4,'k-.', 'linewidth', 2);
xlabel('{\bf problem size n}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
title(sprintf('tic-toc timing averaged over %d runs', nruns),'fontsize',14);
legend('ML, slow evaluation','ML, efficient evaluation',...
    'ML, efficient ev. with reshape','ML, oneliner',...
    'Cpp, slow', 'Cpp fast', 'Cpp superfast', ...
    'O(n^2)','O(n^3)','O(n^4)','location','southeast');
print -depsc2 '../PICTURES/kron_timings_cpp.eps';
