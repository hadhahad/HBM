close all; 
clear all;  %#ok<*CLALL>
syms d m real;
syms om tau al be positive;


%% MODEL %%
% 3rd one: [ (d) <- (m) <- (tau) <- al,be & (d) <- om ]

%%
% Calculation of p(d|m,tau)~N(m,om^(-1):
pd_mtau_nn = exp(-1/2*om*(d-m)^2);   % p(d|m,tau) - not normalized
nc = int(pd_mtau_nn,d,-inf,inf);         % tau is given implicitly
pd_mtau = pd_mtau_nn/nc            % p(d|m,tau) - normalized

%%
% Calculation of p(m|tau)~N(0,tau^(-1)):
pm_tau_nn = exp(-1/2*tau*m^2);    % p(m|tau) - not normalized
nconst = int(pm_tau_nn,m,-inf,inf);      
pm_tau = pm_tau_nn/nc            % p(m|tau) - normalized

%%
% Calculation of p(tau)~Gamma(al,be):
ptau_nn = tau^(al-1)*exp(-be*tau);     % p(tau) - not normalized
nconst_tau = int(ptau_nn,tau,0,inf);
ptau = ptau_nn/nconst_tau            % p(tau) - normalized

%%
% Bayes:
pmtau_d = simplify(pd_mtau*pm_tau*ptau);  % p(m,tau|d) - not normalized
% log:
l_pmtau_d = simplify(log(pmtau_d));     % log(p(m,tau|d)) - not normalized

%%
% Value substitution:
d = 2; al = 1.1; be = 2; om = 1;

%%
% Laplace:
g = gradient(simplify(subs(l_pmtau_d)),[m,tau]);
[m_hat_lap,tau_hat_lap] = vpasolve(subs(g),[m,tau])
m_hat_lap = m_hat_lap(imag(m_hat_lap) == 0);
tau_hat_lap = tau_hat_lap(imag(tau_hat_lap) == 0);

% Hessian evaluation:
H = hessian(subs(l_pmtau_d),[m,tau]);
Sigma = inv(subs(H,[m,tau],[m_hat_lap,tau_hat_lap]));

%%
% Plot Laplace:

vm = d-2:0.01:d+2;                    % range of m
vtau = 0.01:0.001:0.75;                    % range of tau

% figure(1)
% mnorm = normpdf(vm,m_hat_lap,Sigma(1,1));
% plot(vm,mnorm)
% xlabel('m');
% ylabel('p(m|d)');
% title('Laplace, p(m|d)');
% grid

% figure(2)
% taunorm = normpdf(vtau,tau_hat_lap,Sigma(2,2));
% plot(vtau,taunorm)
% xlabel('\tau');
% ylabel('p(\tau|d)');
% title('Laplace, p(\tau|d)');
% grid

%%
% Variational Bayes: 
syms m_hat real;
syms tau_hat var_m al_tau be_tau positive;

l_qm_d = - (m^2*tau)/2 - (om*(d - m)^2)/2;  % ~log(q(m|d))
qm_d = simplify(exp(l_qm_d));              % ~q(m|d)

l_qtau_d = log(tau)*(al - 1)- be*tau - (m^2*tau)/2;      % ~log(q(tau|d))
qtau_d = simplify(exp(l_qtau_d));                            % ~q(tau|d)

% From the above we can conclude, that 
% q(m|d)~N(m_hat_vb, sigma_m) and q(tau|d)~Gamma(al_tau,be_tau). Then:

m_hat_vb = d*om/(tau_hat + om);        % evaluation of m_hat_vb
var_m_vb = 1/(tau_hat + om);                     % evaluation of sigma_m

al_tau_vb = al;                    % evaluation of al_tau
be_tau_vb = (m_hat^2 + var_m)/2 + be;         % evaluation of be_tau

tau_hat_vb = al_tau/be_tau;     % evaluation of tau_hat_vb

%%
% Setting up initial conditions:
tau_hat = 0.5;
var_m = vpa(subs(var_m_vb),5);
m_hat = vpa(subs(m_hat_vb),5);
al_tau = vpa(subs(al_tau_vb));
be_tau = vpa(subs(be_tau_vb),5);

%%
% Iterations:
for i = 1:30
    tau_hat = vpa(subs(tau_hat_vb),5);
    var_m = vpa(subs(var_m_vb),5);
    m_hat = vpa(subs(tau_hat_vb),5);
    be_tau = vpa(subs(be_tau_vb),5);
end

% PROBLEM: Values just start to oscilate...

%%
% Plot Variational Bayes:

vm = d-2:0.01:d+2;                    % range of m
vtau = 0.01:0.001:2;                    % range of tau

figure(3)
mnorm = normpdf(vm,m_hat,var_m);
plot(vm,mnorm)
xlabel('m');
ylabel('p(m|d)');
title('VB, p(m|d)');
grid

figure(4)
taugamma = gampdf(vtau,al_tau,be_tau);
plot(vtau,taugamma)
xlabel('\tau');
ylabel('p(\tau|d)');
title('VB, p(\tau|d)');
grid