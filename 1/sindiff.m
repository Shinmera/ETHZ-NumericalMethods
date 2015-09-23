function sindiff()

x = 1.2;
maxsize = 21;
h = zeros(1,maxsize);
err = zeros(1,maxsize);
msgsize = 0;

for i = 1:maxsize
    % Status message
    fprintf(repmat('\b',1,msgsize));
    msg = sprintf('%3.1f%%\n',i/maxsize*100);
    msgsize = length(msg);
    fprintf('%s',msg);
    % Sample
    h(i) = 10^(1-i);
    f1 = (sin(x+h(i))-sin(x))/h(i);
    f2 = sinderivative(x,h(i));
    err(i) = f1-f2;
end

figure
hold('all')
loglog(h,err,'Color','red');

legend('error')

xlabel('{\bf h}','fontsize',14);
ylabel('{\bf error}','fontsize',14);
print('tex/e2a.png','-dpng');

end