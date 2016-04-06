n = 10;
B = diag(-ones(n-1,1),-1)+diag([2*ones(n-1,1);1],0)...
    + diag(-ones(n-1,1),1);
x = rand(n,1);
fprintf('|x-y| = %d\n',norm(multAmin(B*x)-x));
