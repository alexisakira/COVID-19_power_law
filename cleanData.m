clear
close all
clc;

%% download COVID-19 data and save

% link to data
% COVID-19 cases in US counties from New York Times
url1 = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv';
% Population in US countries from Census Bureau
url2 = 'https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv';

% file name
filename1 = 'us-counties.csv';
filename2 = 'co-est2019-alldata.csv';

% save files using websave function
outfilename1 = websave(filename1,url1);
outfilename2 = websave(filename2,url2);

%% clean COVID-19 data

% get dates
[~,dates] = xlsread('us-counties.csv','A2:A999999'); % import dates
temp = datenum(dates); % convert to datenum
t_unique = unique(temp); % keep unique dates
calendardate = datetime(t_unique,'ConvertFrom','datenum'); % convert to datetime
T = length(calendardate); % number of days

%get counties
fips = xlsread('us-counties.csv','D2:D999999'); % import fips (ID of county)
fips_unique = unique(fips); % keep unique fips
fips_unique = fips_unique(~isnan(fips_unique)); % delete NaN
C = length(fips_unique); % number of counties
[~,county] = xlsread('us-counties.csv','B2:B999999'); % county name

% get cases
cases = xlsread('us-counties.csv','E2:E999999'); % import cases
deaths = xlsread('us-counties.csv','F2:F999999'); % import deaths

covid_cases = zeros(C,T);   % panel of COVID-19 cases
covid_deaths = zeros(C,T);  % panel of COVID-19 deaths

for c=1:C
    for t=1:T
        ind = find((temp==t_unique(t))&(fips==fips_unique(c)));
        if ~isempty(ind)
            covid_cases(c,t) = cases(ind);
            covid_deaths(c,t) = deaths(ind);
        end
    end
end

% adjust for New Yor City
ny_cases = zeros(1,T);
ny_deaths = zeros(1,T);
for t=1:T
    ind = find((temp==t_unique(t))&strcmp(county,'New York City'));
    if ~isempty(ind)
        ny_cases(t) = cases(ind);
        ny_deaths(t) = deaths(ind);
    end
end
covid_cases = [covid_cases;ny_cases];
covid_deaths = [covid_deaths;ny_deaths];

% force cases and deaths to be nondecreasing
covid_cases = cummin(covid_cases,2,'reverse');
covid_deaths = cummin(covid_deaths,2,'reverse');

%% clean population data

state_county = xlsread('co-est2019-alldata.csv','D2:E999999');
pop_estimate = xlsread('co-est2019-alldata.csv','S2:S999999');
fips = 1000*state_county(:,1) + state_county(:,2);

pop = zeros(C,1);
for c=1:C
    ind = find(fips == fips_unique(c));
    pop(c) = pop_estimate(ind);
end
pop = [pop;8.623e6];

C = C+1;
fips_unique = [fips_unique; 99999];

save covid_US_county_data