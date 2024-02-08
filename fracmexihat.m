function [out1,out2, out3] = fracmexihat(LB,UB,N,label)
%MEXIHAT Mexican hat wavelet.
%   [PSI,X] = MEXIHAT(LB,UB,N) returns values of
%   the Mexican hat wavelet on an N point regular
%   grid in the interval [LB,UB].
%   Output arguments are the wavelet function PSI
%   computed on the grid X.
%
%   This wavelet has [-5 5] as effective support.
%
%   See also WAVEINFO.

if nargin < 4,
    label = 'fmxh2.00';
    if nargin < 3,
        N = 168;
        if nargin < 2,
            UB = 8;
            if nargin < 1,
                LB = -5;
            end
        end
    end
end

%[LB,UB,N,label] = deal(varargin{:});
ind = strncmpi('fmxh',label,4);

if isequal(ind,1) label(1:4) = []; end
nu = str2double(deblank(label));

if isempty(nu) error('** fmxh: Invalid wavelet number!'); end
if (nu <= 0) || (nu > 6)
    error('** fmxh: Invalid value for nu **')
end

% Prepare some variables : Mean, StDev and out2
mu      = 10; %(UB - LB)/2;
sigma   = 1; %(UB - LB)/2/5;    % due 5 * sigma
out2    = linspace(LB,UB,N);
t       = out2 + mu;

% Set a default value for nu+v = 2
%if nargin < 5, nu = 2; end

%% Obtain parameters for the calculation of (0 < nu < 1)-order derivative
nux     = nu - floor(nu);
alpha   = nux./(1 - nux);
par1    = -(t - mu - sigma^2.*alpha/2).*alpha;
par3    = (mu/sigma + sigma.*alpha)/sqrt(2);
par2    = t/sigma/sqrt(2) - par3;

% Explicit PseudoGaussian Integrals (mu - x)
fcn     = @(x, a1, t1) exp(-(x.*x - 2*x.*(mu+sigma^2.*a1)+mu^2+2*sigma^2.*a1.*t1)/2/sigma^2);
% for K = 1 : m
Y       = nan(1,N);
for ic = 1 : N,
    Y(1,ic) = integral(@(x) fcn(x,alpha,t(1,ic)), 0, t(1,ic));
end

% Obtain the fractional part
K1      = sqrt(1/2/pi)*(1+alpha);
Df1     = K1.*(exp(par1-par2.^2)-exp(par1-par3.^2))/sigma;
Df2     = K1.*alpha.*Y/sigma;
Df      = Df1-Df2;

% Find first and second ordinary derivatives
dG1      = (exp(-(mu - t).^2/(2*sigma^2)).*(2*mu - 2*t))/(2*sigma^3*sqrt(2*pi));
dG2      = -(exp(-(mu - t).^2/(2*sigma^2)).*(- mu^2 + 2*mu*t + sigma^2 - t.^2))/sigma^5/sqrt(2*pi);
dG3      = -(exp(-(mu - t).^2/(2*sigma^2)).*(mu - t).*(- mu^2 + 2*mu*t + 3*sigma^2 - t.^2))/sigma^7/sqrt(2*pi);
dG4      = (exp(-(mu - t).^2/(2*sigma^2)).*(mu^4 - 4*mu^3*t - 6*mu^2*sigma^2 + 6*mu^2*t.^2 + 12*mu*sigma^2*t - 4*mu*t.^3 + 3*sigma^4 - 6*sigma^2*t.^2 + t.^4))/sigma^9/sqrt(2*pi);
dG5      = (exp(-(mu - t).^2/(2*sigma^2)).*(mu - t).*(mu^4 - 4*mu^3*t - 10*mu^2*sigma^2 + 6*mu^2*t.^2 + 20*mu*sigma^2*t - 4*mu*t.^3 + 15*sigma^4 - 10*sigma^2*t.^2 + t.^4))/sigma^11/sqrt(2*pi);
dG6      = -(exp(-(mu - t).^2/(2*sigma^2)).*(- mu^6 + 6*mu^5*t + 15*mu^4*sigma^2 - 15*mu^4*t.^2 - 60*mu^3*sigma^2*t + 20*mu^3*t.^3 - 45*mu^2*sigma^4 + 90*mu^2*sigma^2*t.^2 - 15*mu^2*t.^4 + 90*mu*sigma^4*t - 60*mu*sigma^2*t.^3 + 6*mu*t.^5 + 15*sigma^6 - 45*sigma^4*t.^2 + 15*sigma^2*t.^4 - t.^6))/sigma^13/sqrt(2*pi);


% Find the high order fractional derivative
switch floor(nu)
    case 1
        dy  = (alpha + 1)*dG1 - alpha*Df;
    case 2
        dy  = (alpha + 1)*((-alpha)*dG1 + dG2) + (- alpha)^2*Df;
    case 3
        dy  = (alpha + 1)*((-alpha)^2*dG1 + (-alpha)*dG2 + dG3) + (- alpha)^3*Df;
    case 4
        dy  = (alpha + 1)*((-alpha)^3*dG1 + (-alpha)^2*dG2 + (-alpha)*dG3 + dG4) + (- alpha)^4*Df;
    case 5
        dy  = (alpha + 1)*((-alpha)^4*dG1 + (-alpha)^3*dG2 + (-alpha)^2*dG3 + (-alpha)*dG4 + dG5) + (- alpha)^5*Df;
    case 6
        dy  = (alpha + 1)*((-alpha)^5*dG1 + (-alpha)^4*dG2 + (-alpha)^3*dG3 + (-alpha)^2*dG4 + (-alpha)*dG5 + dG6) + (- alpha)^6*Df;
end

% Compute values of the Mexican hat wavelet.
out1   = -dy;
out1   = out1/(max(out1));
out3   = dy/max(dy); 
