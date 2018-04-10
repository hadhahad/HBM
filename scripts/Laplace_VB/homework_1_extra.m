close all; clear all %#ok<CLALL>
syms d m real
syms om tau al be positive

%%
pd_mom_nn = exp(-1/2*om*(d-m)^2);   % P(d|m,omega)   nn= not normalized
nk = int(pd_mom_nn,d,-inf,inf);      % int = integrate(function, arg, boundaries)
pd_mom = pd_mom_nn/nk;             % P(d|m,omega) normalized; nk = normalization coefficient

%%
pm_om_nn = exp(-1/2*om*tau*m^2);    % P(m|omega) not normalized
nkm = int(pm_om_nn,m,-inf,inf);     % nkm = normalization coeff.
pm_om = pm_om_nn/nkm;

% the probability distribution of omega is given (a priori) by the gamma
% distribution - we CHOOSE this prior
pom_nn = om^(al-1)*exp(-be*om);     % P(omega) not normalized
nkom = int(pom_nn,om,0,inf);
pom = pom_nn/nkom;

%% Bayes

% simplify = try to make the expression simpler
pmom_d = simplify(pd_mom*pm_om*pom);  % not normalized
lpmom_d = simplify(log(pmom_d));      % not normalized

%% Laplace:                                                          % <------------------------------ HOMEWORK
g = gradient(lpmom_d,[om,m]);
H = hessian(lpmom_d,[om,m]);
[m_hat_lap,om_hat_lap] = solve(g,[m,om])
Sigma = inv(H);

%% Variational Bayes:                                                % <------------------------------ HOMEWORK
syms m_hat real;
syms om_hat var_m al_om be_om positive;

m_vb = d / (1 + tau);
var_vb = 1 / ((1 + tau)*om_hat);
al_vb = al + 1;
be_vb = be + (d^2-2*d*m_hat+(1+tau)*(m_hat^2+var_m))/2;
om_vb = al_om/be_om;


% Initial Conditions:
al = 1; be = 1; tau = 0;     % fixed chosen hyper-parameters
d = 2;                               % measured data point

om_hat = 0.5;
m_hat = subs(m_vb);
var_m = subs(var_vb);
al_om = subs(al_vb);
be_om = subs(be_vb);

% Iterations:
for i = 1:50     % number of iterations was set to 20 for quickness
    om_hat = subs(om_vb);
    m_hat = subs(m_vb);
    var_m = subs(var_vb);
    al_om = subs(al_vb);
    be_om = subs(be_vb);
end

m_hat
om_hat

%% Analytical solution:

tmp_om = simplify(exp(collect(lpmom_d,om))); % temporary omega
tmp_m = simplify(exp(collect(lpmom_d,m)));   % temporary m

% now we can finally integrate
pm_d = int(tmp_om,om,0,inf);          % integrate over omega
pom_d = int(tmp_m,m,-inf,inf);        % integrate over m

lpom_d = simplify(log(pom_d));
tmp_om_d = simplify(exp(collect(lpom_d,m)));   % temporary m
nkom_d = int(tmp_om_d,om,0,inf);

pd = simplify(nkom_d);               % evidence 
% Bayes rule
pm_omd = pmom_d/pom_d;               % conditional posterior

%% Visualize:
al = 1; be = 1; tau = 0;     % fixed chosen hyper-parameters
d = 2;                               % measured data point
vm = d-2:0.1:d+2;                    % range of m
vom = 0.01:0.1:5;                    % range of omega
C = zeros(length(vm),length(vom));   
for i = 1:length(vm)
    m = vm(i);
    for j=1:length(vom)
        om = vom(j);
        C(i,j) = eval(pmom_d/pd);
        C0(i,j) = eval(pm_om*pom);
    end
end

%% Marginal m 
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


%% PLOT:
figure(1)
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
%omnorm = normpdf(vom,eval(om_hat_lap),eval(Sigma(2,2)));
%plot(vom,omnorm,'--')
omgamma = gampdf(vom,al_om,be_om);
plot(vom,omgamma,'--');
hold on
grid
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
%title({'f(\alpha|m)=\int f(x,\alpha) dx'})
legend({'p(\omega|d)','VB'},'Location','northeast')

%side
h3=axes('Position',[0.1, 0.4, 0.2, 0.5]);
hold off
plot(-Em_d,vm);
hold on
%mnorm = normpdf(vm,eval(m_hat_lap),subs(Sigma(1,1)));
%plot(-mnorm,vm,'--')
mnorm = normpdf(vm,m_hat,var_m);
plot(-mnorm,vm,'--');
hold on
grid
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
legend({'p(m|d)','VB'},'Location','north')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0. 7 4];
fig.PaperSize = [7 4];
%print('toy_symb','-dpdf')