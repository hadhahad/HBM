close all; 
clear all;  %#ok<*CLALL>
syms d m real;
syms om tau al be positive;


%% MODEL %%
% 3rd one [ (d) <- (m) <- (tau) <- al,be & (d) <- om ]

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
% Calculation of p(tau)~Gama(al,be):
ptau_nn = tau^(al-1)*exp(-be*tau);     % p(tau) - not normalized
nconst_tau = int(ptau_nn,tau,0,inf);
ptau = ptau_nn/nconst_tau            % p(tau) - normalized

%%
% Bayes:
pmtau_d = simplify(pd_mtau*pm_tau*ptau);  % p(m,tau|d) - not normalized
% log:
l_pmtau_d = log(pmtau_d)     % log(p(m,tau|d)) - not normalized

%%
% Laplace:
d = 2; al = 1.1; be = 1; om = 1;

g = gradient(l_pmtau_d,[m,tau]);
[m_hat,tau_hat] = vpasolve(subs(g),[m,tau])
H = hessian(subs(l_pmtau_d),[m,tau]);
m = m_hat;
tau = tau_hat;
Sigma = (-subs(H))^(-1)
% 
%syms(H,m,m_hat)
%sigma = inv(-H)

%%
% Plot Laplace:

vm = d-2:0.1:d+2;                    % range of m
vom = 0.01:0.1:5;                    % range of omega

figure(1)
mnorm = normpdf(vm,m_hat,Sigma(1,1));
plot(vm,mnorm)
grid

figure(2)
taunorm = normpdf(x,tau_hat,Sigma(2,2));
plot(vm,taunorm)
grid