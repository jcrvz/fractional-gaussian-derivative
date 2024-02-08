function [cfgaussder_normalised,Parameters,t] = cfgaussder_distr (t,nu,mu,sigma)

% Find the start region
new_t               = t(cfgaussder(t,nu,mu,sigma) <= 0);

% Resize t
t                   = linspace(new_t(1),t(end),numel(t));
%mu                  = t(1);

Parameters.Start    = t(1);

% Find the constant
Parameters.Factor   = integral(@(x) cfgaussder(x,nu,mu,sigma).*...
    (-cfgaussder(x,nu,mu,sigma) >= 0), mu,mu + 20*sigma);

% Obtain the distribution
normalised_func         = @(x) cfgaussder(x,nu,mu,sigma).*...
    (-cfgaussder(x,nu,mu,sigma) >= 0)/Parameters.Factor;
cfgaussder_normalised   = normalised_func(t);

% Determine Parameters: Mean, MeanSquare, Variance
Parameters.Norm = integral(@(x) normalised_func(x),...
    mu ,mu + 20*sigma);
Parameters.Mean = integral(@(x) x.*normalised_func(x),...
    mu ,mu + 20*sigma);
Parameters.MeanSquare = integral(@(x) (x.^2).*normalised_func(x),...
    mu ,mu + 20*sigma);
Parameters.Variance = Parameters.MeanSquare - Parameters.Mean^2;
%integral(@(x) ((x - Parameters.Mean).^2).*...
%    normalised_func(x),mu - 20*sigma,mu + 20*sigma);