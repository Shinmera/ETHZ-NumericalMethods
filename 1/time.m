function time()

maxsize = 2;
samples = 10;
n = zeros(1,maxsize);
t1 = zeros(1,maxsize);
t2 = zeros(1,maxsize);
msgsize = 0;

for i = 1:maxsize
    % Status message
    fprintf(repmat('\b',1,msgsize));
    msg = sprintf('%3.1f%%\n',i/maxsize*100);
    msgsize = length(msg);
    fprintf('%s',msg);
    % Sample
    n(i) = i^2;
    for j = 1:samples
        d = rand(n(i),1);
        a = rand(n(i),1);
        x = rand(n(i),1);
        % Calculate and measure
        tic;
        r1 = arrowmatvec(d, a, x);
        t1(i) = t1(i)+toc;
        tic;
        r2 = arrowmatvec2(d, a, x);
        t2(i) = t2(i)+toc;
    end
    t1(i) = t1(i)/samples;
    t2(i) = t2(i)/samples;
end

figure
hold('all')
scatter(n,t1,10,'red','filled');
scatter(n,t2,10,'green','filled');

plot(n,polyval(polyfit(n,t1,2),n),'Color','red');
plot(n,polyval(polyfit(n,t2,2),n),'Color','green');

set(gcf, 'PaperPosition', [0 0 5 3]);
set(gcf, 'PaperSize', [5 3]);
legend('arrowmatvec','arrowmatvec2')

xlabel('{\bf problem size n}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
print('tex/e1e.pdf','-dpdf');

end