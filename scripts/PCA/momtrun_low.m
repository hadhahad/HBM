function [mv,sm,nk]=momtrun_low(m,s,a)
% function evaluates moments of truncated normal distribution
% mv = mean value E(x)
% sm = second moment E(x.^2)
% nk = normalizing coefficient

alf	= (a-m)./(sqrt(2)*s);

% Evaluation has sense for stable only
istab	= find((alf<=3));
astab	= find(alf>3);

%sym mv
mv	= zeros(size(m));
sm	= mv;
nk	= mv;
%%%%%%%% STABLE
if length(istab>0)
	alf	= alf(istab);

	pom	= (1 - erf(alf))* sqrt(pi/2); %delitel
	gam	= ( - exp(-alf.^2)) ./ pom; %phi
	del	= ( - a.*exp(-alf.^2)) ./ pom; %kappa

	mv(istab)	= m(istab) - s(istab).*gam; 
	sm(istab)	= s(istab).^2 + m(istab).*mv(istab) - s(istab).*del;
	nk(istab)	= s(istab) .* pom;
end
%%%%%%%% UNSTABLE
if length(astab)>0
	ma	= m(astab);
	sa	= s(astab);
	mv(astab)	= a - sa.^2./ma;% - (a*sa.^2)./(ma.^2); % - (a^2*sa.^2-2*sa.^4)./(ma.^3);
	sm(astab)	= a^2 - (2*a*sa.^2)./ma + (2*sa.^4 - 2*a^2*sa.^2)./(ma.^2);
	nk(astab)	= -(sa.^2)./ma.*exp(-(a-ma).^2./(2*sa.^2));
end

% check numerics
fault	= find(mv<a);
if ~isempty(fault)
	disp('tN:(mv<a) ');
	keyboard
end
