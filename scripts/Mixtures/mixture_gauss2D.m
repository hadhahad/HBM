function [mu,Sigma]=mixture_gauss2D(d,k,sharp)
% Fitting a mixture of k Gaussians to the input data vector x,
% Prior of each component i the mixture is assumed to be etremely flat:
%   p(mu,sig)= N_mu(0,1/zeta*sig)* G_sig(a0,b0)
% with zeta=0, a0=0, b0=0; 
%
% input:
%   d    - Nx2 vector of data
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

%% ARGUMENTS PROCESSING
if nargin<3
    sharp = 1;
    if nargin < 2
        k = 5;
    end
end

%% INITIALIZATION
data_size = size(d,1);
N = data_size(1);
l = zeros(k,N);

DATA.X = d(:,1);
DATA.Y = d(:,2);

% Split the space between min and max equally
minx = min(d(:,1));
maxx = max(d(:,1));
miny = min(d(:,2));
maxy = max(d(:,2));
distancex = (maxx-minx)/k;
distancey = (maxy-miny)/k;
mu = [(minx + distancex/2 : distancex : maxx); 
    (miny + distancey/2 : distancey : maxy)]' / 10;
Sigma_x = 1/(distancex/sharp)^2;
Sigma_y = 1/(distancey/sharp)^2;
Sigma = diag([Sigma_x, Sigma_y]);
Sigma = Sigma.*ones(2,2,k);

% Assume uniform spread of data between components
priorN = N/k;
alpha = priorN*ones(k,1)/N;

%% ITERATIONS
EM_old = 0;
iter_num = 10000;
for iter = 1:iter_num
    N_distr = zeros(k,1);
    sum_comps = zeros(1,N);
    for comp = 1:k
        l(comp,:) = mvnpdf([DATA.X DATA.Y], mu(comp,:), Sigma(:,:,comp))...
            * alpha(comp);
        sum_comps = sum_comps + l(comp,:);
    end
    l = l./sum_comps;
    N_distr = sum(l,2);
    for comp = 1:k
        sum_mat = zeros(2,2);
        mu(comp,:) = l(comp,:)*[DATA.X DATA.Y]/N_distr(comp);
        for i = 1:N
            sum_mat = sum_mat + ...
                l(comp,i).*(d(i,:)-mu(comp,:))'*(d(i,:)-mu(comp,:));
        end
        Sigma(:,:,comp) = sum_mat./N_distr(comp);
    end
    EM_new = sum(log(sum_comps))
    if abs(EM_new - EM_old) < 1e-3
        break
    end
    EM_old = EM_new;
end

end