%% -- values

% parameters
mu          = 5;
sigma       = sqrt(0.5);

% time series
time_nod    = 500;
time        = linspace(mu,mu + 10*sigma,time_nod);

% nu values
nu_nod      = 11;
n           = 1;
nu          = linspace(0.01 + n - 1,0.9 + n - 1,nu_nod);

% nu            = [0.05:.1:0.9,0.95:0.01:0.99];%[0.05:0.1:0.95];
% nu_nod      = numel(nu);
%% -- calculations
%
% traditional solution
df1dt1      = @(t) (exp(-(mu - t).^2/(2*sigma^2)).*(2*mu - 2*t))/(2*sigma^3*sqrt(2*pi));
df2dt2      = @(t) -(exp(-(mu - t).^2/(2*sigma^2)).*(- mu^2 + 2*mu*t + sigma^2 - t.^2))/sigma^5/sqrt(2*pi);
df3dt3      = @(t) -(exp(-(mu - t).^2/(2*sigma^2)).*(mu - t).*(- mu^2 + 2*mu*t + 3*sigma^2 - t.^2))/sigma^7/sqrt(2*pi);
df4dt4      = @(t) (exp(-(mu - t).^2/(2*sigma^2)).*(mu^4 - 4*mu^3*t - 6*mu^2*sigma^2 + 6*mu^2*t.^2 + 12*mu*sigma^2*t - 4*mu*t.^3 + 3*sigma^4 - 6*sigma^2*t.^2 + t.^4))/sigma^9/sqrt(2*pi);
df5dt5      = @(t) (exp(-(mu - t).^2/(2*sigma^2)).*(mu - t).*(mu^4 - 4*mu^3*t - 10*mu^2*sigma^2 + 6*mu^2*t.^2 + 20*mu*sigma^2*t - 4*mu*t.^3 + 15*sigma^4 - 10*sigma^2*t.^2 + t.^4))/sigma^11/sqrt(2*pi);
df6dt6      = @(t) -(exp(-(mu - t).^2/(2*sigma^2)).*(- mu^6 + 6*mu^5*t + 15*mu^4*sigma^2 - 15*mu^4*t.^2 - 60*mu^3*sigma^2*t + 20*mu^3*t.^3 - 45*mu^2*sigma^4 + 90*mu^2*sigma^2*t.^2 - 15*mu^2*t.^4 + 90*mu*sigma^4*t - 60*mu*sigma^2*t.^3 + 6*mu*t.^5 + 15*sigma^6 - 45*sigma^4*t.^2 + 15*sigma^2*t.^4 - t.^6))/sigma^13/sqrt(2*pi);
%  
% -(t - mu).*exp(-(t - mu).^2/2/sigma^2)/sqrt(2*pi)/sigma^3;

dG1         = df1dt1(time);
dG2         = df2dt2(time);
dG3         = df3dt3(time);
dG4         = df4dt4(time);
dG5         = df5dt5(time);
dG6         = df6dt6(time);
% Gauss-Legendre Approach


% fractional solution
DG          = nan(nu_nod,time_nod);
DGn         = nan(nu_nod,time_nod);

for iinu  = 1 : nu_nod
    nux         = nu(iinu) - floor(nu(iinu));
    gamma_nu    = nux/(1 - nux);
    
    switch floor(nu(iinu))
        case 0
            DG(iinu,:)    = cfgaussder(time,nu(iinu),mu,sigma);
            [DGn(iinu,:),Par(iinu),new_time(iinu,:)] = cfgaussder_distr (time,nu(iinu),mu,sigma);
        case 1
            DG(iinu,:)  = (gamma_nu + 1)*dG1 - gamma_nu*cfgaussder(time,nux,mu,sigma);
        case 2
            DG(iinu,:)  = (gamma_nu + 1)*((-gamma_nu)*dG1 + dG2) + (- gamma_nu)^2*cfgaussder(time,nux,mu,sigma);
        case 3
            DG(iinu,:)  = (gamma_nu + 1)*((-gamma_nu)^2*dG1 + (-gamma_nu)*dG2 + dG3) + (- gamma_nu)^3*cfgaussder(time,nux,mu,sigma);
        case 4
            DG(iinu,:)  = (gamma_nu + 1)*((-gamma_nu)^3*dG1 + (-gamma_nu)^2*dG2 + (-gamma_nu)*dG3 + dG4) + (- gamma_nu)^4*cfgaussder(time,nux,mu,sigma);
        case 5
            DG(iinu,:)  = (gamma_nu + 1)*((-gamma_nu)^4*dG1 + (-gamma_nu)^3*dG2 + (-gamma_nu)^2*dG3 + (-gamma_nu)*dG4 + dG5) + (- gamma_nu)^5*cfgaussder(time,nux,mu,sigma);
    end
end
[T,N]   = meshgrid(time,nu);
%%
new_t       = time(dG1 <= 0);
ord_t       = linspace(new_t(1),time(end),numel(time));
dR          = df1dt1(new_t);



%% -- plot things
%
% parameters
graph.MarkerSize    = 10;
graph.LineWidth     = 1.5;
graph.MinorTicks    = 'on';
graph.FontName      = 'Serif';
graph.FontSize      = 30;

% plotting
imagesize            = [0.9 0.9];
papersize           = [60 20]; % in centimeters

figure('Color','w');%,'Units','Normalized','Position',[0.01 0.05 imagesize]);

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
set(gcf, 'PaperSize', papersize);
set(gcf, 'Position',[papersize.*(1 - imagesize)/2 imagesize.*papersize]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[papersize.*(1 - imagesize)/2 imagesize.*papersize]);
set(gcf, 'Renderer', 'painters');

colours     = [jet(nu_nod);[0 0 0]];
dG_trad     = eval(sprintf('dG%d',n));

%Data = [DG./repmat(max(abs(DG),[],2),1,time_nod);dG_trad./max(abs(dG_trad))]; % normalised
Data  = [DGn;-sqrt(2*pi)*sigma*dR];
stars = [Par.Start,mu];
new_time = [new_time;new_t];

for ii = 1 : nu_nod+1,
    plot(new_time(ii,:)- stars(ii),Data(ii,:),...
        'LineWidth',graph.LineWidth,'Color',colours(ii,:)), %shading interp,
    hold on;
end
box on,

set(gca,'LineWidth',graph.LineWidth,'FontSize',graph.FontSize,'FontName',...
    graph.FontName,'XMinorTick',graph.MinorTicks,'YMinorTick',graph.MinorTicks,...
    'Units','Normalized','Position',[.07 0.1 0.88 0.84],...
    'TickLabelInterpreter', 'latex','YLim',[0 1]),
xlabel('$$\tau = t - \mu_{\nu}$$ (dimensionless)','FontSize',graph.FontSize,'FontName',graph.FontName,...
    'Interpreter','latex');
%ylabel('\nu-order','FontSize',graph.FontSize,'FontName',graph.FontName);
ylabel('$$R_{\sigma}^{\nu}(\tau)$$','FontSize',graph.FontSize,...
    'FontName',graph.FontName,'Interpreter','latex');
lgnd_   = strsplit([sprintf('$$\\nu = %.2f$$,',nu),...
    ['$$\nu = ',sprintf('%.2f',n),'$$,']],',');
lgnd    = lgnd_(1:end-1); xlim([0 5]),
legend(lgnd,'FontSize',graph.FontSize-10,'FontName',graph.FontName,...
    'Interpreter','Latex','Location','Best'),
childrengca = get(gca,'Children');
childrengca(1).LineStyle = '--';

tickValues = get(gca,'YTick');
newLabels = arrayfun(@(value)(sprintf('%.1f',value)), tickValues, 'UniformOutput',false);
set(gca, 'YTickLabel', newLabels);

tickValues = get(gca,'XTick');
newLabels = arrayfun(@(value)(sprintf('%.1f',value)), tickValues, 'UniformOutput',false);
set(gca, 'XTickLabel', newLabels);