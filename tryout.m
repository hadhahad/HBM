close all; 
clear all;  %#ok<*CLALL>
syms d m real;
syms om tau al be positive;


%% MODEL %%
% p(d) ~ N(m, om^(-1))

% Calculation of p(d|m,om)
pd_mtau_nn = exp(-1/2*om*(d-m)^2);   % p(d|m,tau) - not normalized
nc = int(pd_mtau_nn,d,-inf,inf);         % tau is given implicitly
pd_mom = pd_mtau_nn/nc;            % p(d|m,om) - normalized


% Calculation of p(m|tau):
pm_tau_nn = exp(-1/2*om*tau*m^2);    % p(m|tau) - not normalized
nconst = int(pm_tau_nn,m,-inf,inf);      
pm_tau = pm_tau_nn/nc;              % p(m|tau) - normalized


% The probability distribution of 'tau' is set as Gama distribution
% with parameters al, be - determined a priori.
ptau_nn = tau^(al-1)*exp(-be*tau);     % p(tau) - not normalized
nconst_tau = int(ptau_nn,tau,0,inf);
ptau = ptau_nn/nconst_tau;            % p(tau) - normalized

% Bayes:
pm_d = simplify(pd_mtau_nn*pm_tau*ptau);  % p(m,om|d) - not normalized
% log:
l_pmom_d = simplify(log(pm_d));       % log(p(m,om|d)) - not normalized

g = gradient(l_pmom_d,[m,tau]);
H = hessian(l_pmom_d,[m,tau]);
[m_hat,tau_hat] = solve(g,[m,tau])
%syms(H,m,m_hat)
sigma = inv(-H)




%% marginal m 
for i = 1:length(vm)
    m = vm(i);
    om  = 1;
    Em_d(i) = eval(pm_d/pd);
    Em_om_d(i) = eval(pm_omd);
end
% marginal omega
for j=1:length(vom)
    om = vom(j);
    Eom_d(j) = eval(pom_d/pd);
    Eom(j) = eval(pom);
end


%%
figure(2)
% main
h1=axes('Position',[0.4, 0.4, 0.5, 0.5]);
hold off
contour(vom,vm,C,10)
grid
% 	colormap(1-gray);
xlabel('\omega');
ylabel('m');
title(['p(m,om|d=' num2str(d) ')'])
hold on

%below
h2=axes('Position',[0.4, 0.05, 0.5, 0.2]);
hold off
plot(vom,Eom_d);
hold on
plot(vom,Eom,'--');
grid
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
%title({'f(\alpha|m)=\int f(x,\alpha) dx'})
legend({'p(\omega|d)','p(\omega)'},'Location','northeast')
title(' prior vs. posterior')

%side
h3=axes('Position',[0.1, 0.4, 0.2, 0.5]);
hold off
plot(-Em_d,vm);
hold on
plot(-Em_om_d,vm,'--');
hold on
grid
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
legend({'p(m|d)','p(m|om=1,d)'},'Location','north')
title('marginal vs. conditional')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0. 7 4];
fig.PaperSize = [7 4];
%print('toy_symb','-dpdf')