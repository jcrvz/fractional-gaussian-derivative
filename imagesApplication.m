%% -- values
clear all
clc
% parameters
mu          = 10;
sigma       = 1;

% time series
time_nod    = 15;
time        = linspace(max(mu - 3*sigma,0),...
    mu + 3*sigma,time_nod);

% nu values
% nu_nod      = 10;
nu          = [0.5 0.95 1.25 1.5 1.95];
nu_nod      = numel(nu);

%% -- calculations
%
% traditional solution
df0dt0      = @(t) (exp(-(mu - t).^2/(2*sigma^2)))/(sigma*sqrt(2*pi));
df1dt1      = @(t) (exp(-(mu - t).^2/(2*sigma^2)).*(2*mu - 2*t))/(2*sigma^3*sqrt(2*pi));

dG0         = df0dt0(time);
dG1         = df1dt1(time);

% Components
dG0x        = repmat(dG0,time_nod,1);

% fractional solution
DG          = nan(nu_nod,time_nod);
for iinu  = 1 : nu_nod
    nux         = nu(iinu) - floor(nu(iinu));
    gamma_nu    = nux/(1 - nux);
    
    switch floor(nu(iinu))
        case 0
            DG(iinu,:)    = cfgaussder(time,nu(iinu),mu,sigma);
        case 1
            DG(iinu,:)  = (gamma_nu + 1)*dG1 - gamma_nu*cfgaussder(time,nux,mu,sigma);
   end
end
[T,N]   = meshgrid(time,nu);

%% Image Applications
Im  = imread('flor.jpg');
Im2 = single(rgb2gray(Im));
for iinu = 1 : nu_nod
    figure('Name',sprintf('nu = %.4f',nu(iinu))),
    
    switch floor(nu(iinu))
        case 0
            Filtro = double(DG(iinu,:));
            Im3 = conv2(Im2, Filtro ,'same');
            Im4 = conv2(Im2, Filtro','same');
            Im5 = sqrt((Im3).^2+(Im4).^2);
        case 1
            DG_x = repmat(DG(iinu,:),time_nod,1);
            Filtro = DG_x.*dG0x' + DG_x'.*dG0x;
            Im5 = imfilter(Im2, Filtro);%,'same');
            Im5 = Im5./sum(Im5(:));
    end
    
    MMax = max(max(Im5));
    MMin = min(min(Im5));
    filtered_image = uint8(255*(Im5-MMin)/(MMax-MMin));
    
    J = imresize(filtered_image, 0.15);
    imwrite(J,sprintf('filtered_image%d.png',iinu));
    % imwrite(imresize(uint8(Im2),0.15),sprintf('filtered_image%d.png',0));
    imshow(J);
end