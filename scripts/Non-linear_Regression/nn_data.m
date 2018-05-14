x = -10:0.1:10;
y0= 1./abs(x).*sin(abs(x));

ind = sort(ceil(rand(1,50)*length(x)));
x0 = x(ind);
y = y0(ind) + 0.01*randn(size(ind));

save nn_data x0 y y0 

fig = figure;
plot(x0,y,'x')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0. 4 3];
fig.PaperSize = [4 3];
print('nn_data','-dpdf')