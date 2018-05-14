function [mu,sig,L,Smu,DEBUG]=mixture_gauss(x,k,sharp)
% Fitting a mixture of k Gaussians to the input data vector x,
% Prior of each component i the mixture is assumed to be etremely flat:
%   p(mu,sig)= N_mu(0,1/zeta*sig)* G_sig(a0,b0)
% with zeta=0, a0=0, b0=0; 
%
% input:
%   x    - Nx1 vector of data
%   k    - maximum number of components
%   sharp- initial sharpness of the components, relative to the spread of
%          the data. The variance of initial components is:
%             std(x)_k = (max(x)-min(x))/k/sharp
%          Higher values of "sharp" favors smaller clusters.
%
% output:
%   mu   - mean values of the components
%   sig  - variances of the components
%   L    - probabilities of membership in each component
%   Smu  - variances of the mean value estimates (var(mu))

if nargin<3
    sharp = 1
    if nargin<2
        k = 5;
    end
end

N = size(x,1);
L = zeros(N,k);

%% INITIALIZATION
minx = min(x);
maxx = max(x);

X = x*ones(1,k);
ON = ones(N,1);

% split the space between min and max equally
distance = (maxx-minx)/k;
mu = (minx+distance/2: distance : maxx)'/10;
om = 1/(distance/sharp)^2 * ones(k,1);

% assume uniform spread of data between components
priorN = N/k;
alpha = priorN*ones(k,1);
beta = alpha./om;
Smu = om/priorN;
a0= 1e-10;
b0= 1e-10;
w = ones(1,k)/k;

%% Iterations
Lold = L;
Mu = ON*mu';

%% DEBUG
niter =900;
DEBUG.mu = zeros(k,niter);
DEBUG.sig = zeros(k,niter);
DEBUG.sL = zeros(k,niter);

for ite = 1:niter

    % new L
    Om = ON*om';
    Lam = (-(X-Mu).^2.*Om + ON*log(om)' ); % log p()
    Lnn = exp(0.5*Lam-max(max(Lam)))*diag(w); % exp(log P())* w

    
    L = Lnn ./(sum(Lnn,2)*ones(1,k)+1e-10);
    sumL = sum(L,1)'+1e-10;
    
    % new mu
    mu = (L'*x)./(sumL+1e-10);
    Mu = ON*mu';

    % new Om
    alpha = a0 + sumL +1;
    beta = b0 + sum(L.*((X-Mu).^2),1)';
    
    om = alpha./beta;
    Smu = 1./(sumL.*om);
    
    % new weight
    w = sumL/sum(sumL);
    
    %te
    if 0
        figure(1);
        subplot(2,1,1);
        [hy,hx]=hist(x,100);

        sL = sum(L,1);
        fig=figure(1);
        hold off
        bar(hx,hy);
        hold on
        for kk=1:k
            fest = exp(-0.5*(hx-mu(kk)).^2*om(kk));
            f = w(kk)*fest/sum(fest)*N;
            plot(hx,f,'LineWidth',3);
        end
        title(['EM algorithm iteration ' num2str(ite)])
        
        subplot(2,1,2);
        hold off
        stem(x,L(:,1))
        hold on
        stem(x,L(:,2),'r')
        fig.PaperPosition = [-0.3 0 4.3 3];
        fig.PaperSize = [4 3];
        %         print(['mix_em' num2str(ite) ],'-dpdf')
        % pause
    end

    
    if 1
        DEBUG.mu(:,ite) = mu;
        DEBUG.sig(:,ite) = 1./om;
        DEBUG.sL(:,ite) = sumL;

    end
    if norm(Lold-L)<1e-3
        break
    end


    Lold=L;
end

%% output
sig = 1./om;
DEBUG.mu=DEBUG.mu(:,1:ite);
DEBUG.sig=DEBUG.sig(:,1:ite);
DEBUG.sL = DEBUG.sL(:,1:ite);
