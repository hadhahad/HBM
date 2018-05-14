clc;
clear all; %#ok<CLALL>

m= [0,4,4];
sig = [1,1,1];
N = [200,100,100];

x = [m(1)+sqrt(sig(1))*randn(N(1),1);
     m(2)+sqrt(sig(2))*randn(N(2),1);
     m(3)+sqrt(sig(3))*randn(N(3),1)];

%%
K=3;
[mu,sig,L,Smu,DBG]=mixture_gauss1D(x,K,4)


[hy,hx]=hist(x,40);

sL = sum(L,1);
figure(1); 
hold off
bar(hx,hy);
hold on
for k=1:K
    fest = exp(-0.5*(hx-mu(k)).^2/sig(k));
    f = sL(k)*fest/sum(fest);
    plot(hx,f,'LineWidth',3);
end