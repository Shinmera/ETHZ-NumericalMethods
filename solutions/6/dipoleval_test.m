function dipoleval_test()

%%% testing script for dipoleval
n = 10;  % number of interpolation points
N = 100; % number of evaluation points
t = linspace(0,1,n);
y = rand(1,n);
x = linspace(0,1,N);

dp = dipoleval(t,y,x);
P = polyfit(t,y,n-1);
p = polyval(P,x);
dp_alt = dipoleval_alt(t,y,x);

figure;
hold on;
plot(t,y,'ro');
plot(x,p,'b-');
plot(x,dp,'g-.');
plot(x,dp_alt,'r--');
legend('data points','polynomial','AN-derivative','naive-derivative',...
       'location','best');
print -depsc '../PICTURES/dipoleval_test.eps'
hold off;

end

function dp = dipoleval(t,y,x)
dp=zeros(1,length(x));

for i = 1:length(x)
    n  = length(y);
    p  = y;
    dP = zeros(1,n);
    for im=2:n
        for i0=im-1:-1:1
            % compute dp(i)'s
            dP(i0) = p(i0+1) + (x(i)-t(i0))*dP(i0+1) - p(:,i0) - (x(i)-t(im))*dP(i0);
            dP(i0) = dP(i0) / ( t(im) - t(i0) );
            % compute p(i)'s
            p(i0)  = (x(i)-t(i0))*p(i0+1) - (x(i)-t(im))*p(i0);
            p(i0)  = p(i0) / ( t(im) - t(i0) );
        end
    end
    dp(i)=  dP(1);
end

end

function dp = dipoleval_alt(t,y,x)

n  = length(y);
P = polyfit(t,y,n-1);
dP = linspace(n-1,0,n).*P;
dp = polyval(dP(1:end-1),x);

end
