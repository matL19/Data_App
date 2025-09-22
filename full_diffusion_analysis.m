%% FULL DATA ANALYSIS of diffusion experiment - from start to finish -

% ---- Variables that require user manual input will be surrounded ----
% like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.

% ----
date_of_experiment = "2025-09-18";
% ----

year_of_experiment = year(datetime(date_of_experiment));

% ---- here put your path to CHEM-SGR (end with a slash). This will be necessary
% throughout the script.
isilon_path = "/Volumes/CHEM-SGR/";
% ----
data_path = isilon_path + "sgr-ftir.chem.pitt.edu/" + year_of_experiment + "/" + date_of_experiment;

try
    cd(isilon_path)
catch
    warning("An error ocurred accessing the data directory. Make sure you are connected to Isilon.")
end
%% Step 1: process the microscope image

%% Load the image from Isilon

% ---- path to the image within CHEM-SGR (end with a slash) ----
pre_image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-09-18/";
% ----

% ---- name of the pre-diffusion image file (use the .tif one) ----
pre_image_filename = "PMNTF2 EMIM FTIR windows 20250918 65 C pre-diffusion.tif";
% ----

%% Process the image

% Get the radius from the image
% ------------
dx_def = "from center";
rad_def = "circle fit";
units = "mm";
path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
ref_filename = "scaled_reticle01.tif";
% ------------

[pre_radius,pre_displacement,pre_other] = radius_from_image(isilon_path+pre_image_path+pre_image_filename,...
    path_to_ref,ref_filename,...
    "image scale",80,"dx definition",dx_def,"radius definition",rad_def,...
    "units",units);

fprintf("r = %4f mm; dx = %4f mm\n",pre_radius,pre_displacement)

%% Process the post image

% ---- is the run complete and you have a post image?
post_image_complete = true;
% ----

if post_image_complete
    
    % ---- path to the image within CHEM-SGR (end with a slash) ----
    post_image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-09-19/";
    % ----
    
    % ---- name of the pre-diffusion image file (use the .tif one) ----
    post_image_filename = "PMNTF2 EMIM FTIR windows 20250918 65 C post-diffusion.tif";
    % ----
    
    % Get the radius from the image
    % ------------
    [post_radius,post_displacement,post_other] = radius_from_image(isilon_path+post_image_path+post_image_filename,...
        path_to_ref,ref_filename,...
        "image scale",40,"dx definition",dx_def,"radius definition",rad_def,...
        "units",units);
    % ------------

    
end

%% Show the difference in pre and post image
if post_image_complete
    % Show the image comparison
    pre_image = imread(isilon_path + pre_image_path + pre_image_filename);
    post_image = imread(isilon_path + post_image_path + post_image_filename);
    [optimizer, metric] = imregconfig('multimodal');
    post_image_moved = imregister(post_image,pre_image,'similarity',optimizer,metric);
    figure(5015);
    imshowpair(pre_image,post_image_moved)
end
%% Save the data
% make a neat structure
image_processing.pre.path = pre_image_path + pre_image_filename;
image_processing.pre.radius = pre_radius;
image_processing.pre.displacement = pre_displacement;
image_processing.pre.units = units;
image_processing.pre.dx_definition = dx_def;
image_processing.pre.radius_definition = rad_def;
image_processing.pre.other_data = pre_other;

if post_image_complete
    image_processing.post.path = post_image_path + post_image_filename;
    image_processing.post.radius = post_radius;
    image_processing.post.displacement = post_displacement;
    image_processing.post.units = units;
    image_processing.post.dx_definition = dx_def;
    image_processing.post.radius_definition = rad_def;
    image_processing.post.other_data = post_other;
end

cd(data_path)
save(date_of_experiment + "_image_processing_data.mat",'image_processing')
%% ---- END OF STEP 1 ----

%% Step 2: Fit the FTIR peaks to obtain the uptake curve

%% Load in the spectra
cd ~
% --- the indicies of the spectra you wish to use ----
spectra_range = [1:316];
% ----

% data to use
% -------
use_spectra = [1:316];
% -------
% some of the data can be bad. only these spectra will be used for all data
% analysis and visualization.

% --- the spectra file prefix ---
file_prefix = 'PMNTF2EMIM_FTIR_20250918_65C_';
% ----

% --- the name of the temperature file
temperature_log_filename = "TemperatureLog[12_30_29_PM][9_18_2025].log";
% ---

% --- experimental parameters ---
volume = NaN;  % in microliters
path_length = 67.5083;  % in microns
gel_radius = pre_radius*1000;  % in microns
displacement = pre_displacement*1000; % in microns
time_delay = 300;  % between spectra, in seconds
sample_name = "PMIM NTF2";
your_name = "Matt";
temperature_setpoint = 65;
% ---

cd(data_path)
[data1,freq] = LoadSpectra(data_path,file_prefix,spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end

% Subtract to the initial spectrum
sub_data = data1 - data1(:,1);

% INITIALIZE OBJECT
f = FTIRexperiment(sub_data,freq,volume,path_length,gel_radius,...
    time_delay,sample_name,date_of_experiment,your_name);
f = f.timeAxis(data_path,file_prefix,spectra_range);
f.displacement = displacement;
f.temperature_setpoint = temperature_setpoint;

% Determine the temperature from the log file
temp_log = readmatrix(temperature_log_filename);
if ~isempty(temp_log)
    temp_log = temp_log(:,2);
    temp_log = temp_log(~isnan(temp_log));
    f.temperature = mean(temp_log);
    f.temperature_std = std(temp_log);
    fprintf("Temperature recorded: %.3f ± %.3f ºC\n",f.temperature,f.temperature_std)
else
    fprintf("Cannot read temperature log file until experiment is complete.\n")
end

fprintf("Successfully imported " + size(f.data,2) + " spectra.\n")

%% Guesses for FTIR peak fitting, by eye
% ---- Which spectrum will you match to? Usually the last one is good.
trial_spectrum = 316;
% ----

% set the fit range. Usually doesn't need to be changed
range1 = [2290 2390];

clear sp
% ---- User-input starting point values ----
sp.center = 2340;
sp.wg = 1.7;
sp.wl = 1.7;
sp.a1 = 3.87;  % main peak height
sp.a2 = 0.07; % expected Boltzmann factor for bend
sp.a3 = 0; % gas lines
sp.c0 = 0.008;
sp.c1 = 0; % baseline slope
% ----

%fit function requires fliipped inputs
freq = flip(f.freqAxis);
s = flip(f.data(:,trial_spectrum));


%get x and y data for the fit
ind1 = find(freq>=range1(1) & freq<range1(2));
x = freq(ind1);
ydata = s(ind1);

%plot the fitted function using user parameters
yfit = co2GasLineFitFunction(x,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1);
res = ydata-yfit;
sse = sum(res.^2);

figure(182);clf
plot(x,ydata,'o',x,yfit,x,res-0.1,'r-o')
%% Do the FTIR peak fit
T = tic; %time the fitting for later display
f = gasLineFit(f,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1,...
    "TolFun",1e-8,"TolX",1e-8,"MaxFunEvals",5000);
stop = toc(T);

%selecte 4 evenly placed fits to plot
n_spectra = size(f.data,2);
iis = ceil([1 n_spectra/4 n_spectra*3/4 n_spectra]);
figure(209);clf
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
%% Plotting the uptake curve for viewing
figure(319);clf

for ii = use_spectra
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

time_axis = f.timePts;
concovertime = f.concOverTime;
plot(subplot(2,1,2),time_axis(use_spectra),concovertime(use_spectra),'o-','color','blue');
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

%% Water content analysis

% wavenumber range of the water peaks
% -------
region_of_interest = [3150 3900];
% -------

% integration range
% -------
integration_range = [3300 3650];
% -------

% frequency point to track
% -------
freqpt = 3513;
% -------

% Plot of water over time

y_max = zeros(1,numel(use_spectra));
y_int = y_max;
t = f.timePts/3600; % convert to hours

figure(4745);clf
clear sub_freq sub_data
subplot(2,1,2)
hold on

for ii = 1:numel(t)
    
    [sub_data(:,ii),sub_freq] = getDataSubset(f.freqAxis,f.data(:,ii),region_of_interest);
    
    sub_data(:,ii) = baselineCorrect(sub_freq,sub_data(:,ii),"value",3200);
    
    curve_color = [ii/numel(t) 0 (numel(t)-ii)/numel(t)];
    
    if any(use_spectra == ii)
        plot(sub_freq,sub_data(:,ii),'Color',curve_color)
    end
    
    y_max(ii) = getDataSubset(sub_freq,sub_data(:,ii),[freqpt freqpt]);
    
    y_int(ii) = trapz(getDataSubset(sub_freq,sub_data(:,ii),integration_range));
    
end
freqline = xline(freqpt,'g--','LineWidth',2);
int_range = plot([integration_range(1) integration_range(2)],[max(sub_data(:,end)) max(sub_data(:,end))],'go--','LineWidth',2);
xlim(region_of_interest)
xlabel('Frequency (cm^{-1})')
ylabel('Absorbance')
legend([freqline, int_range],{'Frequency point to track','Integration range'})
box off
hold off

subplot(2,1,1)
yyaxis left
plot(t(use_spectra),y_max(use_spectra),'bo-','MarkerFaceColor','blue','MarkerEdgeColor','blue')
ylabel('Max Absorbance (O.D.)')
yyaxis right
plot(t(use_spectra),y_int(use_spectra),'ro-','MarkerFaceColor','red','MarkerEdgeColor','red')
ylabel('Integrated Intensity (A.U.)')
xlabel('Time (hr)')
legend('Absorbance at frequency point','Integrated intensity','Location','northwest')
box off

set(gcf,'Position',[946     1   648   946])
set(gcf,'Color','white')

%% Show second derivative
figure(7331);clf
hold on
plot(f.timePts,gradient(gradient(f.concOverTime)),'b-')
plot(f.timePts,movmean(gradient(gradient(f.concOverTime)),10),'ro-')
plot([0 max(f.timePts)],[0 0],'r--')
legend('d^2c/dt^2','10-pt moving avg.','Location','northeast')
title('Second derivative of conc over time plot')
xlabel('Time (s)')
ylabel('d^2 Concentration / dt^2')
box off
set(gca,'TickDir','out')
hold off
%% Save the data
FTIR_peak_fitting_params.sp = sp;
FTIR_peak_fitting_params.trial_spectrum = trial_spectrum;
FTIR_peak_fitting_params.fit_range = range1;
FTIR_peak_fitting_params.spectra_range = spectra_range;
FTIR_peak_fitting_params.file_prefix = file_prefix;
cd(data_path)
save(date_of_experiment + "_FTIR_peak_fitting.mat",'FTIR_peak_fitting_params')
%% ---- END OF STEP 2 ----

%% Step 3: Fit for single diffusion coefficient

%% Guesses for uptake curve fitting, by eye
t = f.timePts;
t = t(use_spectra);
y = f.concOverTime;
y = y(use_spectra);
A = f.radius;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dy = 0;

% ---- User input starting values
dx = f.displacement;  % from the image analysis. can be overridden
%     D      C
sp = [22.48  0.07225]; % put guess here
ub = [1e5 10];
lb = [0 0];
% ----

figure(728);clf
plot(t,y)
hold on
ymodel = diffusion_moving_beam(t,sp(1),A,sp(2),nmax,sigma,dx,dy,"rlim",rlim,'dx definition','from center');
plot(t,ymodel)
res = y(:) - ymodel(:);
plot(t,res-0.025,'ro')
errs(1) = sum((res).^2);
errs(2) = std(res);
fprintf("SSE: " + errs(1) + "; Std. res: " + errs(2) + "\n")
%% Do the uptake curve fitting

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter',...
    'Weights',ones(size(y))*1/errs(2),...
    'TolFun',1e-16,...
    'TolX',1e-16);

ft = fittype(@(D,C,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,"rlim",rlim,'dx definition','from center'),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D','C'},...
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
figure(144);clf

plot(f.diffusionFitResult.x,f.diffusionFitResult.ydata,...
    'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(f.diffusionFitResult.x,f.diffusionFitResult.yfit,...
    'red','LineWidth',1.5)
residuals = f.diffusionFitResult.ydata(:) - f.diffusionFitResult.yfit(:);
plot(f.diffusionFitResult.x,(residuals*1 - 0.02),'o','MarkerEdgeColor','red')
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
%% Save the data

single_diffusion_fitting_params.sp = sp;
single_diffusion_fitting_params.ub = ub;
single_diffusion_fitting_params.lb = lb;
single_diffusion_fitting_params.dx = dx;
single_diffusion_fitting_params.nmax = nmax;
single_diffusion_fitting_params.rlim = rlim;
single_diffusion_fitting_params.sigma = sigma;
save(date_of_experiment + "_single_diffusion_fitting_params.mat",'single_diffusion_fitting_params')

f.fitMethod = 'diffusion_moving_beam.m';
cd(data_path);
save(date_of_experiment + "_single_diffusion_fitting_params.mat",'single_diffusion_fitting_params')
save(f.dateString + "_single_diffusion_fitting","f")

%% Update lab notebook with results

% ---- PUT IN THE CORRECT NOTEBOOK PAGE TITLE ---
notebook = 'Matt Lab Notebook';
folder = 'Experiments';
page_title = '2025-09-18 Diffusion of CO2 in PMNTF2 EMIM at 65 C';
% ----

obj = labarchivesCallObj('notebook',notebook,...
    'folder',folder,...
    'page',page_title);

% pre photo
% ----------
figure(1)
% ----------
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + pre_radius + "mm, " + "dx = " + pre_displacement + "mm.";
obj = obj.updateFigureAttachment('caption',caption);

if post_image_complete
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + post_radius + "mm, " + "dx = " + post_displacement + "mm.";
% post photo
% ----------
figure(2)
% ----------
obj = obj.updateFigureAttachment('caption',caption);

end

% comparison photo
% ----------
figure(5015)
% ----------
obj = obj.updateFigureAttachment('caption',"Comparison photo between the pre- and post-diffusion image. Differences are highlighted in green and pink. Images were translated and rotated to find best match using imregconfig.");

% uptake curve
% ----------
figure(319)
% ----------
obj = obj.updateFigureAttachment;

% single diffusion fitting result
% ----------
figure(144)
% ----------
caption = "Single diffusion coefficient fitting: ";
coeffs = coeffnames(f.diffusionFitResult.fobj);
units = ["um^2/s" "M"];
if numel(units) ~= numel(coeffs)
    error("Cannot match all fitting parameters with a unit.")
end
for ii = 1:numel(coeffs)
    std_devs{ii} = (ci(2,ii) - ci(1,ii))/4;
    if isnan(std_devs{ii})
        caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
            + " " + units(ii) + ", ";
    else
        caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
            + " ± " + std_devs{ii} + " " + units(ii) + ", ";
    end
end
obj = obj.updateFigureAttachment('caption',caption);

% water content
% ----------
figure(4745)
% ----------
obj = obj.updateFigureAttachment('caption',"Water content throughout the experiment.");

%% ---- END OF STEP 3 ----

%% Step 4: Fit for double diffusion coefficient

% Proceed to the double_diffusion_fitting.m script for this part. You will
% need to have run this entire script for it to work as it will draw from
% the .mat files you saved here.