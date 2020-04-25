%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randgrowth_covid.m
% (c) 2020 Brendan K. Beare and Alexis Akira Toda
% 
% Purpose: 
%       Replication file of the analysis in "On the Emergence of a Power
%       Law in the Distribution of COVID-19 Cases" by Brendan K. Beare and
%       Alexis Akira Toda
%
% Version 1.0: April 25, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cleanData % run this once to clean data
clear
clc;
load covid_US_county_data % load data file (unnecessary if run cleanData.m)

%% figure formatting

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')

set(0,'DefaultTextFontSize', 14)
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultLineLineWidth',2)

colors = get(gca,'ColorOrder');

c1 = colors(1,:);
c2 = colors(2,:);
close all

%% set estimation period
t_0 = datetime([2020,01,21]); % start date of data
t_end = datetime([2020,03,31]); % up to March 2020
T_incl = days(t_end-t_0) + 1; % number of days to include

%% define variables
growth = log(covid_cases(:,2:T_incl)) - log(covid_cases(:,1:T_incl-1)); % log growth rate
T_epidemic = sum(covid_cases(:,1:T_incl)>0,2); % days since epidemic

%% Kolmogorov-Smirnov test of identical growth distribution between two dates (not in paper)

S = T_incl-2; % number of tests
Pval_KS = nan(1,S);
for s=1:S
    x1 = growth(:,s);
    x1 = x1(x1>0);
    x1 = x1(x1<Inf);

    x2 = growth(:,s+1);
    x2 = x2(x2>0);
    x2 = x2(x2<Inf);
    if (~isempty(x1))&&(~isempty(x2))
        [~,p] = kstest2(x1,x2);
        Pval_KS(s) = p;
    end
end

figure
plot(calendardate(1:S),Pval_KS);
xlabel('Date')
ylabel('P-value of Kolmogorov-Smirnov test')

sum(~isnan(Pval_KS))
sum(Pval_KS<0.05)

%% OLS estimation to test random growth
Coeff = nan(4,S);
SE = nan(4,S);
Ssize = zeros(1,S);
N_min = 30; % set minimum number of observations in cross-section
ind_incl = []; % index of dates to include
for s=1:S
    ind = find((growth(:,s)>=0)&(growth(:,s)<Inf)&(growth(:,s+1)>=0)&(growth(:,s+1)<Inf)...
        &(covid_cases(:,s)>0)); % index of valid data points
    if length(ind)>= N_min
        ind_incl = [ind_incl s];
        y = growth(ind,s+1); % log growth rate
        temp = sum(covid_cases(:,1:s)>0,2); % days since epidemic
        X = [growth(ind,s) log(covid_cases(ind,s)) temp(ind)]; % regressors
        mdl = fitlm(X,y); % OLS estimation
        Coeff(:,s) = mdl.Coefficients.Estimate; % coefficients
        SE(:,s) = mdl.Coefficients.SE; % standard errors
        Ssize(s) = length(ind); % sample size
    end
end

titleList = {'Constant','Growth','Size','Days'};

% plot results

%set(0,'DefaultTextFontSize', 22)
%set(0,'DefaultAxesFontSize', 22)

figure
plot(calendardate(ind_incl+1),Coeff(1,ind_incl),'-','Color',c1); hold on
plot(calendardate(ind_incl+1),Coeff(1,ind_incl)+1.96*SE(1,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),Coeff(1,ind_incl)-1.96*SE(1,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),zeros(1,length(ind_incl)),'k-','LineWidth',1); hold off
ylabel('Estimate')
legend('Point estimate','95\% confidence band')
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
title(titleList{1})
%{
%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_OLS_' titleList{1}],'-dpdf')
%}

figure
plot(calendardate(ind_incl+1),Coeff(2,ind_incl),'-','Color',c1); hold on
plot(calendardate(ind_incl+1),Coeff(2,ind_incl)+1.96*SE(2,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),Coeff(2,ind_incl)-1.96*SE(2,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),zeros(2,length(ind_incl)),'k-','LineWidth',1); hold off
xlabel('Date')
ylabel('Estimate')
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
title(titleList{2})
%{
%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_OLS_' titleList{2}],'-dpdf')
%}

figure
plot(calendardate(ind_incl+1),Coeff(3,ind_incl),'-','Color',c1); hold on
plot(calendardate(ind_incl+1),Coeff(3,ind_incl)+1.96*SE(3,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),Coeff(3,ind_incl)-1.96*SE(3,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),zeros(3,length(ind_incl)),'k-','LineWidth',1); hold off
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
title(titleList{3})
%{
%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_OLS_' titleList{3}],'-dpdf')
%}

figure
plot(calendardate(ind_incl+1),Coeff(4,ind_incl),'-','Color',c1); hold on
plot(calendardate(ind_incl+1),Coeff(4,ind_incl)+1.96*SE(4,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),Coeff(4,ind_incl)-1.96*SE(4,ind_incl),':','Color',c1);
plot(calendardate(ind_incl+1),zeros(4,length(ind_incl)),'k-','LineWidth',1); hold off
xlabel('Date')
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
title(titleList{4})
%{
%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_OLS_' titleList{4}],'-dpdf')
%}

%set(0,'DefaultTextFontSize', 14)
%set(0,'DefaultAxesFontSize', 14)


%% analysis of age distribution

% discrete logistic distribution

trunc = 1; % indicator of truncation
Tmax = max(T_epidemic);
counts = zeros(1,Tmax+1);
for t=0:Tmax
    counts(t+1) = sum(T_epidemic == t); % age distribution
end

tic
[param,PMF] = discreteLogistic_ML(counts(2:end),trunc);
toc
p_event = param(2); % birth probability

if trunc == 0
    temp = T_epidemic;
else
    temp = T_epidemic(T_epidemic>0);
end

figure
histogram(temp,'normalization','pdf'); hold on
plot([trunc:Tmax],PMF); hold off
xlabel('Days since first confirmed case')
ylabel('Probability mass')
legend('Data','Truncated logistic fit')
%{
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_hist_Days','-dpdf')
%}

%% analysis of growth distribution
Nmin = 10; % lower bound to include
gvec = growth(:);
temp = covid_cases(:,1:T_incl-1);
cvec = temp(:);
ind1 = find((~isnan(gvec)&(gvec > 0))&(gvec<Inf)&(cvec >= Nmin));
% include only positive growth for nicer picture
ind0 = find((~isnan(gvec)&(gvec >= 0))&(gvec<Inf)&(cvec >= Nmin));
% include zero growth
N0 = length(ind0) - length(ind1); % number of zero growth observations
pi = N0/length(ind0); % fraction of zero growth

[phat,pci] = gamfit(gvec(ind1)); % fit gamma distribution
alpha = phat(1);
beta = 1/phat(2);

temp = (1-p_event)*(1-pi)/(1-pi+p_event*pi);
zeta = beta*(1-temp^(1/alpha)); % theoretical Pareto exponent

% plot determination of Pareto exponent
nGrid = 100;
zGrid = linspace(0,beta,nGrid+1);
zGrid(end) = [];
Mz = pi + (1-pi)*(1-zGrid/beta).^(-alpha);
q = 1-p_event;

z1 = 1.5;
M1 = pi + (1-pi)*(1-z1/beta).^(-alpha);
delta = 0.02;

figure
plot(zGrid,Mz); hold on
plot(zeta*ones(1,2),[0 1/q],'k--','LineWidth',1);
text(zeta+delta,delta,'$\zeta$','HorizontalAlignment','left','VerticalAlignment','bottom');
plot(zGrid,1/q*ones(1,nGrid),'k--','LineWidth',1);
text(delta,1/q+delta,'$1/q$','HorizontalAlignment','left','VerticalAlignment','bottom');
text(z1,M1,'$M(z)$','HorizontalAlignment','right','VerticalAlignment','bottom');
title('Determination of Pareto exponent')
ylim([0 2])
xlim([0 2])
xlabel('$z$')

% plot histogram of growth rate distribution
x = linspace(0,max(gvec(ind1)),1000);

figure
histogram(gvec(ind1),'normalization','pdf'); hold on
plot(x,gampdf(x,phat(1),phat(2)),'Color',c2); hold off
xlabel('Growth rate of confirmed cases')
ylabel('Probability density')
legend('Data','Gamma fit')
%{
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_hist_Growth','-dpdf')
%}

%% Pareto exponent
% cases
case_end = covid_cases(:,T_incl); % number of cases at end of sample
case0 = zeros(C,1); % number of cases at beginning of sample
for c=1:C
    ind = find(covid_cases(c,:)>0,1,'first');
    if ~isempty(ind)
        case0(c) = covid_cases(c,ind);
    end
end

%temp = case_end./case0; % growth in cases

temp = case_end; % delete this and replace with previous line if desired

% inspect largest 10 counties
size_fips = [temp fips_unique];
ind = find(~isnan(size_fips(:,1)));
size_fips = size_fips(ind,:);
[~,ind] = sort(size_fips(:,1),'descend');
size_fips = size_fips(ind,:);
size_fips(1:10,:) % largest 10 counties

temp = sort(temp(temp>0),'descend');

N = length(temp); % sample size

[alpha,xmin] = plfit(temp); % use Clauset-Shalizi-Newman code
zetahat = alpha-1; % estimated Pareto exponent - CSN uses unusual convention
N_tail = sum(temp>=xmin);
se_zetahat = zetahat/sqrt(N_tail); % standard error from Hill estimator
%Pval = plpva(temp,xmin); % P-value of Kolmogorov-Smirnov test of Pareto distribution

% log-log plot
figure
loglog(temp,[1:N]/N,'o'); hold on;
loglog(temp(1:N_tail),(N_tail/N)*(temp(1:N_tail)/temp(N_tail)).^(-zetahat),'-','Color',c2);
xlabel('Confirmed COVID-19 cases on March 31')
ylabel('Tail probability')
legend('Data','Pareto fit')
xticks([10^0,10^1,10^2,10^3,10^4,10^5]);
text(temp(floor(N_tail/2)),1/10,['$\hat{\zeta}=$' num2str(zetahat,'%4.3f')...
    ' (' num2str(se_zetahat,'%4.3f') ')'],'HorizontalAlignment','left');
text(temp(1)-5000,1/(N+1000),['$\nearrow$'],'HorizontalAlignment','right');
text(temp(1)-22500,1/(N+1750),['New York City'],'HorizontalAlignment','right');
text(temp(2)-1000,2/(N+1000),['$\nearrow$'],'HorizontalAlignment','right');
text(temp(2)-5000,2/(N+1750),['Nassau County'],'HorizontalAlignment','right');
text(temp(6)-300,6/(N+1000),['$\nearrow$'],'HorizontalAlignment','right');
text(temp(6)-1800,6/(N+1750),['Los Angeles County'],'HorizontalAlignment','right');
%{
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_Pareto_Cases','-dpdf')
%}

% population (not in paper)
temp = sort(pop,'descend');
temp = temp(temp > 1e5);
N = length(temp);

[alpha,xmin] = plfit(temp);
zetahatpop = alpha-1;
N_tail = sum(temp>=xmin);
se_zetahatpop = zetahatpop/sqrt(N_tail); % standard error

figure
loglog(temp,[1:N]/N,'o'); hold on;
loglog(temp(1:N_tail),(N_tail/N)*(temp(1:N_tail)/temp(N_tail)).^(-zetahatpop),'-','Color',c2);
xlabel('Population')
ylabel('Tail probability')
legend('Data','Pareto fit')
text(temp(floor(N_tail/2)),1/10,['$\widehat{\zeta}=$' num2str(round(zetahatpop,3))...
    ' (' num2str(round(se_zetahatpop,3)) ')'],'HorizontalAlignment','left');
%{
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_Pareto_Population','-dpdf')
%}