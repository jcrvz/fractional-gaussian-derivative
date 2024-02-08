function Df = cfgaussder (t,nu,mu,sigma)

% Read input arguments
%
if nargin < 4
    sigma = sqrt(2)/2;
    if nargin < 3
        mu = 0;
        if nargin < 2
            nu = 0.5;
            if nargin < 1
                error('There is not value for t, Df = cfgaussder(t)');
            end
        end
    end
end

if size(t,1) ~= 1, t = t'; end

%% CRUZ and AVINA :: MODIFIED
% Calculate the theoretical function value
alpha   = nu./(1 - nu);
par1    = -(t - mu - sigma^2.*alpha/2).*alpha;
par3    = (mu/sigma + sigma.*alpha)/sqrt(2);
par2    = t/sigma/sqrt(2) - par3;

%% Explicit PseudoGaussian Integrals (mu - x).
[rows,columns] = size(t);
fcn = @(x, a1, t1) exp(-(x.*x - 2*x.*(mu + sigma^2.*a1) + mu^2 + ...
    2*sigma^2.*a1.*t1)/2/sigma^2);

Y   = nan(rows,columns);
for ir = 1 : rows
    for ic = 1 : columns
        Y(ir,ic) = integral(@(x) fcn(x,alpha,t(ir,ic)), 0, t(ir,ic));
    end
end

K1 = sqrt(1/2/pi)*(1+alpha);
Df1 = K1.*(exp(par1-par2.^2)-exp(par1-par3.^2))/sigma;
Df2 = K1.*alpha.*Y/sigma;
Df  = Df1-Df2;
