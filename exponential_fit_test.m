
plot(x, exp(-x/0.1))

x=[1,2,4,6,8]';
y=[100,140,160,170,175].';
g = fittype('a-b*exp(-x/tau)');
% g = fittype('a-b*exp(-c*x)');
% f0 = fit(x,y,g,'StartPoint',[-100 -5])
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1])
% f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(1,8,50);

figure;
plot(f0,'r-');

figure;
plot(x,y,'o',xx,f0(xx),'r-');

figure;
hold on;
plot(x,y,'o')
plot(xx,f0(xx),'r-');
hold off;

figure;
hold on;
plot(x,y,'o')
plot(f0,'r-');
hold off;

figure;
hold on;
plot(x,y,'o')
plot(x,f0(x),'r-');
hold off;

figure;
hold on;
plot(x,y,'o')
plot(x,f0,'r-');
hold off;


% ones(size(x))
% -exp(-x)
% [ones(size(x)), -exp(-x)]\y
% [[ones(size(x)), -exp(-x)]\y; 1]

