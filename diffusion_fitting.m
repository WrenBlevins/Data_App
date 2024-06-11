%% This script is for MANUAL FITTING and DATA ANALYSIS
% of the FTIR data for use OUTSIDE of the app.
%%
[data1,freq] = LoadSpectra("/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2024/2024-06-11",...
    "ML_20240611_",[2:38]);
%% set up object
f = FTIRexperiment(data1,freq(:,1),0.2,12,3140,300,"70% EMIM NTF2 in PEGDA","2024-06-11","Matt");
% f = f.timeAxis;
f = f.timeAxis("/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2024/2024-06-11",...
    "ML_20240611_",[2:38]);
%% make initial guesses
% have the user select which spectrum to guess from
ii = 37;

% set the fit range
range1 = [2290 2390];

% set starting point using values from the user
center = 2339;
wg = 1.9;
wl = 1.7;
a1 = 0.55;  % main peak height
a2 = 0.07; % expected Boltzmann factor for bend
a3 = 0.008; % gas lines
c0 = 0.8;
c1 = 1e-4; % baseline slope

%fit function requires fliipped inputs
freq = flip(f.freqAxis);
s = flip(f.data(:,ii));


%get x and y data for the fit
ind1 = find(freq>=range1(1) & freq<range1(2));
x = freq(ind1);
ydata = s(ind1);

%plot the fitted function using user parameters
yfit = co2GasLineFitFunction(x,center,wg,wl,a1,a2,a3,c0,c1);
res = ydata-yfit;
sse = sum(res.^2);

figure(1);clf
plot(x,ydata,'o',x,yfit,x,res-0.1,'r-o')
%app.UIAxes3.Title = (sprintf('Initial guess SSE = %f',sse));
%% do the gas line fit
T = tic; %time the fitting for later display
f = gasLineFit(f,center,wg,wl,...
    a1,a2,a3,c0,...
    c1);
stop = toc(T);

%selecte 4 evenly placed fits to plot
n_spectra = size(f.data,2);
iis = ceil([1 n_spectra/4 n_spectra*3/4 n_spectra]);
figure(2);clf
for ii = iis
    plot(f.fittedSpectra(ii).x,f.fittedSpectra(ii).ydata,'o',...
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).yfit,...
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).res-0.1,'ro')
    hold on
end
hold off

%let the user know how it went
review = "";
tl = 0;
for ii = 1:n_spectra
    if f.fittedSpectra(ii).O.exitflag < 1
        review = [review;'Spectrum '+ii+' did not converge!!! Results might not be trustworthy.'];
        tl = tl+1;
    end
end
if tl==0
    review = "All fits were successful.";
end
review = [review;"Fitting took "+stop+" seconds."];
review
%% plotting the fits
figure(3);clf

% number of spectra to show
n = size(f.data,2);

%find the indicies for the amount of spectra desired
spectraIndicies = zeros(1,n);
interval = ceil(size(f.data,2)/n);
for ii = 1:n
    spectraIndicies(ii) = (ii*interval);
end

for ii = spectraIndicies
    temp = f.fittedSpectra(ii).fobj;
    pf = co2GasLineFitFunction(f.fittedSpectra(ii).x,temp.center,temp.w_g,temp.w_l,...
        temp.a1,temp.a2,0,0,0);
    plot(subplot(2,1,1),f.fittedSpectra(ii).x,pf)
    hold on
end
title('Fitted Spectra')
xlabel('Wavenumbers (cm^{-1})')
ylabel('Absorbance (AU)')
box off
set(gca,'TickDir','out')
hold off

%converts data in f to vector of concentration values of CO2
%for each spectrum, generates time axis, returns a plot of
%concentration vs time
%syntax: plotConcOverTime(f)
plot(subplot(2,1,2),f.timePts,concOverTime(f),'color','blue');
hold on
title('Concentration Over Time')
xlabel('Time (s)')
ylabel('Concentration (M)')
box off
set(gca,'TickDir','out')
hold off

set(gcf,'Units','normalized')
set(gcf,'Color','w')
set(gcf,'Position',[0.5 0 0.35 1])
%% final conc if applicable

% f.finalSpectrum = LoadSpectra();
% f = f.getFinalConc;

%% fit for diffusion coefficient
%get parameters ready
t = f.timePts;
%         t = t(1:end-3);
%         t = t-t(1);
y = f.concOverTime;
%         y = y(4:end);
A = f.radius;
C = f.finalConc;
nmax = 150;
rres = 50;
rlim = 700;
sigma = 704;
dy = 0;
sp = [150 0.220 0]; % put guess here
ub = [1e5 1 0.5*f.radius];
lb = [0 0 0];

figure(728);clf
plot(t,y)
hold on
plot(t,diffusion_moving_beam(t,sp(1),f.radius,sp(2),nmax,sigma,sp(3),dy))


%%

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter');

ft = fittype(@(D,C,dx,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D','C','dx'},...
    'options',opts);

%set up structure for storing output
out = struct('x',[],'ydata',[],'yfit',[],'res',[],...
    'fobj',[],'G',[],'O',[]);

tic

%do the fit
[fobj,G,O] = fit(t,y',ft);

toc
%%
%get results
yfit = fobj(t);
out.x = t;
out.ydata = y;
out.yfit = yfit;
out.res = y - yfit;
out.fobj = fobj;
out.G = G;
out.O = O;

if out.O.exitflag < 1
    warning('Curve fit did not converge!!! Results might not be trustworthy.');
end

figure(4);clf

plot(out.x,out.ydata,'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(out.x,out.yfit,'red','LineWidth',1.5)
residuals = out.yfit - out.ydata(:);
plot(out.x,(residuals*10 - 0.02),'o','MarkerEdgeColor','red')
legend('Data points','Fitted curve','Location','northwest')
hold off


% get confidence intervals
ci = confint(out.fobj);

readout = [string(out.fobj.D)]
others = ["95% Confidence Interval is "+ci(1)+" to "+ci(2)+".";...
    "R^2 = "+string(out.G.rsquare)]

fobj
%%
figure(712);clf
plot(g.timePts,g.concOverTime,'blue')
hold on
plot(f.timePts,f.concOverTime,'red')
% plot(f.timePts,f.concOverTime,'black')
legend("today","the other day",'Location','northwest')