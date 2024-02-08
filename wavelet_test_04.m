%% Create (renew) pattern data

% renew the patterns' folder
namedir = 'patterns-fgd';
eval(['!rm -rv ', namedir]);
eval(['!mkdir ', namedir]);
results_path = './ecg100m';

% n is the total derivative order
for n = 1.125 : 0.125 : 5.75
    namef = sprintf('fmxh%g',n*1000);
    named = sprintf('fmxh%.3f',n);
    [inverse_pattern,~,forward_pattern] = fracmexihat(-5,8,168,named);
    save(['patterns-fgd/',namef,'.mat'],...
        'inverse_pattern','forward_pattern');
    fprintf('File %s.mat created!\n',namef);
    %plot(pattern), getframe(gcf);
end

fontSize = 14;
%% Load experimental ECG

% plotATM('100m');

% Load ECG data from MIT-100m at physionet.org/cgi-bin/atm/ATM
[time,signal_MLII,signal_V5] = plotATM('100m',false);

% Make the signal shorter and good-to-compute
%lpow = floor(log2(length(time)));
samples = 1 : numel(time);
s = signal_V5(samples);
t = time(samples);

s = (s - min(s))/diff(minmax(s));

% Set the malfunction
t_initial = 6.5*0.2*4 + 0.2;
t_final = t_initial + 0.2*4;

fig = figure('color','white','name',...
    sprintf('Signal V5'),'Unit','Normalized',...
    'Position',[0.1 0.05 0.8 0.8],'MenuBar','none',...
    'PaperOrientation','landscape','PaperUnit','inch','PaperSize',[11 8.5/3]);
hold on;
fill([t_initial t_initial t_initial t_final t_final t_final],[0 0 1 1 0 0],[1 .85 .85],'LineStyle','none');
h1 = plot(t,s,'b','linewidth',1.5);
xlim([0 10]); ylim([0 1.0]);
%title(sprintf('V5 signal MIT-100m'),...
%    'Interpreter','LaTeX','FontSize',20);
xlabel('Time [s]','Interpreter','LaTeX','FontSize',fontSize,'LineWidth',1.5);
ylabel('Norm. Amplitude',...
    'Interpreter' ,'LaTeX','FontSize',fontSize,'LineWidth',1.5);
legend({'Reported Arrhythmia','Signal from V5'},'Interpreter','LaTeX','FontSize',fontSize)
set(h1.Parent,'TickLabelInterpreter','LaTeX','FontSize',fontSize,...
    'LineWidth',1.5,'Box','on');
    set(gca,'xticklabel',num2str(num2str(get(gca,'xtick')','%.0f')));
    %set(gca,'yticklabel',num2str(num2str(get(gca,'ytick')','%.2f')));
h1.Parent.XAxis.TickDirection = 'out';
h1.Parent.YAxis.TickDirection = 'out';

%print(fig,[results_path,'/','ECG100m_','signal'],'-dpdf','-r300','-fillpage');
%    fprintf(sprintf('Figure %s generated!\n',['signal']));
%close(fig);

%% Perform the computations (SWT)

% Prepare folders and files
%!mkdir ecg100m
files = dir([namedir,'/']);
namefiles_ = {files.name};
namefiles  = namefiles_(3:end);
nnamefiles = numel(namefiles);

% Number of decomposition levels
levels = 3;

% Load patterns from the folder
DATAa = nan(nnamefiles,numel(time));
DATAd = nan(nnamefiles,numel(time));
orders = nan(1,nnamefiles);
orderslabel = cell(1,nnamefiles);
for ii = 1 : nnamefiles
    filestr = namefiles{ii};
    load([namedir,'/',filestr]);
    X = linspace(0,1,numel(forward_pattern));
    Y = forward_pattern;
    %     Y = inverse_pattern;
    
    [psi,xval,nc] = pat2cwav(Y, 'orthconst',3,'continuous') ;
    
    % Split the file name in two strings
    headstr = filestr(1:4);
    tailstr = filestr(5:8);
    
    % Get the fractional order
    fracorder = str2double(tailstr)/1000;
    
    % Get the continous wavelet transform
    %[cfs,f] = cwt(s,psi,Fs);
    
    % Obtain filters
    [Lo_D,Hi_D,Lo_R,Hi_R] = orthfilt(psi);
    
    % Calculate swt
    [swa,swd] = swt(s,levels,Lo_D,Hi_D);
    
    % Store these values in a matrix
    DATAa(ii,:) = (swa(levels,:)-min(swa(levels,:)))/diff(minmax(swa(levels,:)));%(swa(levels,:)-min(swa(levels,:)))/diff(minmax(swa(levels,:)));
    %DATAd(ii,:) = swd(levels,:);%/max(abs(swd(levels,:)));
    orders(ii) = fracorder;
    orderslabel = ['$$\\nu+n=$$ ',sprintf('%.2f',fracorder)];
    
    % Plot it in pretty good resolution
    fig = figure('color','white','name',...
        sprintf('Frac MEX Hat, nu = %.2f',fracorder),'Unit','Normalized',...
        'Position',[0.1 0.05 0.8 0.4],'MenuBar','none',...
        'PaperOrientation','landscape','PaperUnit','inch','PaperSize',[12 8.5/3]);
        
    y_min = min(DATAa(ii,:)); y_max = max(DATAa(ii,:)); hold on;
    fill([t_initial t_initial t_initial t_final t_final t_final],[y_min 0 y_max y_max 0 y_min],[1 .85 .85],'LineStyle','none');
    plot(t,DATAa(ii,:),'b','linewidth',1.5);
    xlim([0 10]); ylim([y_min y_max]);
    %title(sprintf('Level %d coefficients',levels),...
    %    'Interpreter','LaTeX','FontSize',20);
    xlabel('Time [s]','Interpreter','LaTeX','FontSize',fontSize,'LineWidth',1.5);
    ylabel(sprintf('Norm. Coeff. $$\\hat{a}_%d $$',levels),...
        'Interpreter','LaTeX','FontSize',fontSize,'LineWidth',1.5);
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',fontSize,...
        'LineWidth',1.5,'Box','on');
    set(gca,'xticklabel',num2str(num2str(get(gca,'xtick')','%.0f')));
    %set(gca,'yticklabel',num2str(num2str(get(gca,'ytick')','%.2f')));
    h1 = get(gca); 
    h1.XAxis.TickDirection = 'out';
    h1.YAxis.TickDirection = 'out';
    
    print(fig,[results_path,'/','ECG100m_',headstr,tailstr],'-dpdf','-r300','-fillpage');
    fprintf(sprintf('Figure %s generated!\n',[headstr,tailstr]));    
    close(fig);
end

%% Summary plot

fig = figure('color','white','name',...
        sprintf('Frac MEX Hat, nu = %.2f',fracorder),'Unit','Normalized',...
        'Position',[0.1 0.05 0.8 0.8],'MenuBar','none',...
        'PaperOrientation','landscape','PaperUnit','inch','PaperSize',[12  8.5/2]);
        
    h0 = boxplot(DATAa',orders,'OutlierSize',10,...
        'Symbol','r.','Widths',0.3,'Whisker',1);
    xlabel('($$n+\nu$$)','Interpreter','LaTeX','FontSize',fontSize,'LineWidth',1.5);
    ylabel(sprintf('Norm. Coeff. $$\\hat{a}_%d $$',levels),...
        'Interpreter','LaTeX','FontSize',fontSize,'LineWidth',1.5);
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',fontSize,...
        'LineWidth',1.5,'Box','on');
    xtickangle(90)
    h1 = get(gca); 
    h1.XAxis.TickDirection = 'out';
    h1.YAxis.TickDirection = 'out';
    set(gca,'xticklabel',num2str(num2str(orders','%.3f')));
    %set(gca,'yticklabel',num2str(num2str(get(gca,'ytick')','%.2f')));
    set(findobj(gca,'type','line'),'linew',1.5)
    
    print(fig,[results_path,'/','ECG100m_','summary'],'-dpdf','-r300','-fillpage');
    fprintf(sprintf('Figure %s generated!\n','summary'));    
    close(fig);