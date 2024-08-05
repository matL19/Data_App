%% This script is for MANUAL FITTING and DATA ANALYSIS
% of the FTIR data for use OUTSIDE of the app.
%%
cd ~  % This is here because sometimes MATLAB gets confused 
% finding the Isilon folder so you have to reset the current folder to
% somewhere on the disk first.
spectra_range = [2:77];
[data1,freq] = LoadSpectra(...
    '/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2024/2024-08-02',...
    "ML_20240802_2_",spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end
% [data1,freq] = LoadSpectra();

% INITIALIZE OBJECT
f = FTIRexperiment(data1,freq,0.08,12,1550,120,"Pure EMIM NTf2","2024-08-02","Matt");
% f = f.timeAxis;
f = f.timeAxis(...
    '/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2024/2024-08-02',...
    "ML_20240802_2_",spectra_range);
% f = f.timeAxis;

clear spectra_range
%% make initial guesses
% have the user select which spectrum to guess from
ii = 60;

% set the fit range
range1 = [2290 2390];

% set starting point using values from the user
center = 2342.5;
wg = 1.7; 
wl = 1.7;
a1 = 2.7;  % main peak height
a2 = 0.07; % expected Boltzmann factor for bend
a3 = 0.02; % gas lines
c0 = 1.18;
c1 = 5.5e-5; % baseline slope

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
plot(x,ydata,'o',x,yfit,x,res+1,'r-o')
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
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).res + 1,'ro')
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

plot(subplot(2,1,2),f.timePts,concOverTime(f),'o-','color','blue');
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
%           t = t(1:end-15);
%         t = t-t(1);
y = f.concOverTime;
%         y = y(4:end);
%           y = y(1:end-15);
A = f.radius;
C = f.finalConc;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dx = 0;
dy = 0;
sp = [200 0.3 0]; % put guess here
ub = [1e5 1e3 0.5*f.radius];
lb = [0 0 0];

figure(728);clf
plot(t,y)
hold on
plot(t,diffusion_moving_beam(t,sp(1),f.radius,sp(2),nmax,sigma,sp(3),dy,"rlim",rlim))


%% Actually do the fit

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter');

ft = fittype(@(D,C,dx,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,"rlim",rlim),...
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

f.diffusionFitResult = out;
%% display fit result
figure(4);clf

plot(f.diffusionFitResult.x,f.diffusionFitResult.ydata,...
    'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(f.diffusionFitResult.x,f.diffusionFitResult.yfit,...
    'red','LineWidth',1.5)
residuals = f.diffusionFitResult.yfit - f.diffusionFitResult.ydata(:);
plot(f.diffusionFitResult.x,(residuals*10 - 0.02),'o','MarkerEdgeColor','red')
legend('Data points','Fitted curve','Location','northwest')
xlabel('Time (s)')
ylabel('Concentration (M)')
hold off


% get confidence intervals
ci = confint(f.diffusionFitResult.fobj);

readout = [string(f.diffusionFitResult.fobj.D)]
others = ["95% Confidence Interval is "+ci(1)+" to "+ci(2)+".";...
    "R^2 = "+string(f.diffusionFitResult.G.rsquare)]

f.diffusionFitResult.fobj

%% Update lab notebook with results
f.fitMethod = 'diffusion_moving_beam.m';
cd("/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2024/"+f.dateString);
save(f.dateString,"f")

obj = labarchivesCallObj('notebook','Matt Lab Notebook',...
    'folder','Experiments',...
    'page','2024-08-02 Measurement of Diffusion Coefficient of CO2 in pure EMIM NTf2');
figure(3)
obj = obj.updateFigureAttachment;
figure(4)
obj = obj.updateFigureAttachment('caption',...
    "D = "+fobj.D+" um^2/s,C = "+fobj.C+" M,dx = "+fobj.dx+" um");

%%
this_file_name = 'diffusion_fitting.m';
new_file_name = 'your_date_string_here.m';
copyfile("diffusion_fitting.m","/Users/matthewliberatore/Desktop/")
