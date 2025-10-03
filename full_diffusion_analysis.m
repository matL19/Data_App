%% FULL DATA ANALYSIS of diffusion experiment - from start to finish -

% ---- Variables that require user manual input will be surrounded ----
% like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.

% ----
date_of_experiment = "2025-08-27";
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
%%

% ========================================================
%   Step 1: Process the image
% ========================================================


%% Load the image from Isilon

% ---- path to the image within CHEM-SGR (end with a slash) ----
pre_image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-08-27/";
% ----

% ---- name of the pre-diffusion image file (use the .tif one) ----
pre_image_filename = "PMIM NTF2 2D windows 20250827 rt pre-diffusion.tif";
% ----

%% Process the image

% Get the radius from the image
% ------------
dx_def = "from edge";
rad_def = "circle fit";
units = "mm";
path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
ref_filename = "scaled_reticle01.tif";
% ------------

[pre_radius,pre_displacement,pre_other] = radius_from_image(isilon_path+pre_image_path+pre_image_filename,...
    path_to_ref,ref_filename,...
    "image scale",80,"dx definition",dx_def,"radius definition",rad_def,...
    "units",units,"flag_plot",false);

fprintf("r = %4f mm; dx = %4f mm\n",pre_radius,pre_displacement)

% make a neat structure
image_processing.pre.path = pre_image_path + pre_image_filename;
image_processing.pre.radius = pre_radius;
image_processing.pre.displacement = pre_displacement;
image_processing.pre.units = units;
image_processing.pre.dx_definition = dx_def;
image_processing.pre.radius_definition = rad_def;
image_processing.pre.other_data = pre_other;

%% Save the data - just pre image

cd(data_path)
save(date_of_experiment + "_image_processing_data.mat",'image_processing')

%% Show the resulting pre-diffusion image analysis

if ~exist('image_processing','var')
    cd(data_path)
    load(date_of_experiment + "_image_processing_data.mat")
end
I = imread(isilon_path + image_processing.pre.path);
figure(1);clf
plotGelImageAnalysis(I*100,image_processing.pre.radius,image_processing.pre.displacement,image_processing.pre.other_data,image_processing.pre.units)

%% Process the post image

% ---- is the run complete and you have a post image?
post_image_complete = true;
% ----

if post_image_complete
    
    % ---- path to the image within CHEM-SGR (end with a slash) ----
    post_image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-10-02/";
    % ----
    
    % ---- name of the pre-diffusion image file (use the .tif one) ----
    post_image_filename = "PMIM NTF2 2D windows 20251001 80 C post-diffusion.tif";
    % ----
    
    path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
    ref_filename = "scaled_reticle01.tif";
    
    % Get the radius from the image
    % ------------
    [post_radius,post_displacement,post_other] = radius_from_image(isilon_path+post_image_path+post_image_filename,...
        path_to_ref,ref_filename,...
        "image scale",40,"dx definition",image_processing.pre.dx_definition,...
        "radius definition",image_processing.pre.radius_definition,...
        "units",image_processing.pre.other_data.units,"flag_plot",false);
    % ------------
    
    image_processing.post.path = post_image_path + post_image_filename;
    image_processing.post.radius = post_radius;
    image_processing.post.displacement = post_displacement;
    image_processing.post.units = image_processing.pre.other_data.units;
    image_processing.post.dx_definition = image_processing.pre.dx_definition;
    image_processing.post.radius_definition = image_processing.pre.radius_definition;
    image_processing.post.other_data = post_other;
    
end

%% Save the data again - with post image

if post_image_complete
    
    cd(data_path)
    save(date_of_experiment + "_image_processing_data.mat",'image_processing')
    
end

%% Show the post image analysis

if post_image_complete
    
    if ~exist('image_processing','var')
        cd(data_path)
        load(date_of_experiment + "_image_processing_data.mat")
    end
    I = imread(isilon_path + image_processing.post.path);
    figure(2);clf
    plotGelImageAnalysis(I*100,image_processing.post.radius,image_processing.post.displacement,image_processing.post.other_data,image_processing.post.units)
    
end

%% Show the difference in pre and post image

if post_image_complete
    
    % Show the image comparison
    pre_image = imread(isilon_path + image_processing.pre.path);
    post_image = imread(isilon_path + image_processing.post.path);
    [optimizer, metric] = imregconfig('multimodal');
    post_image_moved = imregister(post_image,pre_image,'similarity',optimizer,metric);
    figure(5015);
    imshowpair(pre_image,post_image_moved)
    
end
%% ---- END OF STEP 1 ----

%%

% ========================================================
%   Step 2: Fit the FTIR peaks to obtain the uptake curve
% ========================================================


%% Load in the spectra
cd ~
% --- the indicies of the spectra you wish to use ----
spectra_range = [1:483];
% ----

% --- the spectra file prefix ---
file_prefix = 'PMIMNTF2_2D_20251001_80C_';
% ----

% --- the name of the temperature file
temperature_log_filename = "TemperatureLog[3_13_50_PM][10_1_2025].log";
% ---

% --- experimental parameters ---
volume = NaN;  % in microliters
path_length = 10.7889;  % in microns
time_delay = 120;  % between spectra, in seconds
sample_name = "PMIM NTF2";
your_name = "Matt";
temperature_setpoint = 80;
% ---

% load in image processing data if we do not have it
if ~exist('image_processing','var')
    cd(data_path)
    load(date_of_experiment + "_image_processing_data.mat")
end
% set the displacement and radius
if isfield('image_processing','post')
    gel_radius = mean([image_processing.pre.radius image_processing.post.radius]);
    displacement = mean([image_processing.pre.displacement image_processing.post.displacement]);
else
    gel_radius = image_processing.pre.radius;
    displacement = image_processing.pre.displacement;
end

% convert to microns
if image_processing.pre.units == "mm"
    gel_radius = gel_radius * 1000;
    displacement = displacement * 1000;
end

% load in the data
cd(data_path)
tic
[data1,freq] = LoadSpectra(data_path,file_prefix,spectra_range);
T = toc;
fprintf("Time to load %i spectra: %.3f mins\n",size(data1,2),T/60);
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
try
    temp_log = readmatrix(temperature_log_filename);
catch
    warning("No temperature log found. A temperature was not recorded.")
    temp_log = [];
end
if ~isempty(temp_log)
    temp_log = temp_log(:,2);
    temp_log = temp_log(~isnan(temp_log));
    f.temperature = mean(temp_log);
    f.temperature_std = std(temp_log);
    fprintf("Temperature recorded: %.3f ± %.3f ºC\n",f.temperature,f.temperature_std)
else
    fprintf("Cannot read temperature log file until experiment is complete.\n")
end


fprintf("Successfully imported %i spectra.\n",size(f.data,2))


%% Guesses for FTIR peak fitting, by eye
% ---- Which spectrum will you match to? Usually the last one is good.
trial_spectrum = 483;
% ----

% set the fit range. Usually doesn't need to be changed
range1 = [2290 2390];

clear sp
% ---- User-input starting point values ----
sp.center = 2342.25;
sp.wg = 1.7;
sp.wl = 1.7;
sp.a1 = 0.6;  % main peak height
sp.a2 = 0.07; % expected Boltzmann factor for bend
sp.a3 = 0.01; % gas lines
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

%% Package up the data and save it

% peak fitting
FTIR_peak_fitting_params.sp = sp;
FTIR_peak_fitting_params.trial_spectrum = trial_spectrum;
FTIR_peak_fitting_params.fit_range = range1;
FTIR_peak_fitting_params.spectra_range = spectra_range;
FTIR_peak_fitting_params.file_prefix = file_prefix;
cd(data_path)
save(date_of_experiment + "_FTIR_peak_fitting.mat",'FTIR_peak_fitting_params')

% FTIR object
cd(data_path)
save(date_of_experiment + "_single_diffusion_fitting.mat","f")

%% Plotting the uptake curve for viewing

if ~exist('f','var')
    cd(data_path)
    load(date_of_experiment + "_single_diffusion_fitting.mat")
end

% data to use
% -------
use_spectra = [1:size(f.data,2)];
% -------
% some of the data can be bad. only these spectra will be used for all data
% analysis and visualization.


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

%% ---- END OF STEP 2 ----

%%

% ========================================================
%       Step 3: Fit for single diffusion coefficient
% ========================================================

%% Make sure we have data

% FTIR object
if ~exist('f','var')
    cd(data_path)
    load(date_of_experiment + "_single_diffusion_fitting.mat")
end

% image analysis
if ~exist('image_processing','var')
    cd(data_path)
    load(date_of_experiment + "_image_processing_data.mat")
end

%% Guesses for uptake curve fitting, by eye

% ---- User input starting values
dx = f.displacement;
%     D      C
sp = [65  0.073]; % put guess here
ub = [1e5 10];
lb = [0 0];
use_spectra = [1:size(f.data,2)];
% ----

% time axis
t = f.timePts;
t = t(use_spectra);

% concentration data
y = f.concOverTime;
y = y(use_spectra);

% other model parameters
A = f.radius;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dy = 0;
dx_def = image_processing.pre.dx_definition;

figure(728);clf
plot(t,y)
hold on
ymodel = diffusion_moving_beam(t,sp(1),A,sp(2),nmax,sigma,dx,dy,"rlim",rlim,'dx definition',dx_def);
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

ft = fittype(@(D,C,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,"rlim",rlim,'dx definition',dx_def),...
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

 %% Save the data

% diffusion fitting parameters
single_diffusion_fitting_params.sp = sp;
single_diffusion_fitting_params.ub = ub;
single_diffusion_fitting_params.lb = lb;
single_diffusion_fitting_params.dx = dx;
single_diffusion_fitting_params.nmax = nmax;
single_diffusion_fitting_params.rlim = rlim;
single_diffusion_fitting_params.sigma = sigma;
single_diffusion_fitting_params.dx_def = dx_def;

cd(data_path)
save(f.dateString + "_single_diffusion_fitting_params.mat",'single_diffusion_fitting_params')

f.fitMethod = 'diffusion_moving_beam.m';
cd(data_path);
save(f.dateString + "_single_diffusion_fitting","f")

%% display fit result

if ~exist('f','var')
    cd(data_path)
    load(date_of_experiment + "_single_diffusion_fitting.mat")
end

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

% show cfit object
f.diffusionFitResult.fobj

% display values and standard deviations
ci = confint(f.diffusionFitResult.fobj);
D_std = (ci(2,1) - ci(1,1))/4;
C_std = (ci(2,2) - ci(1,2))/4;
fprintf("\n============\nResults\n============\n\tD = %.4f ± %.4f μm^2/s\n\tC = %.4f ± %.4f M\n\n",f.diffusionFitResult.fobj.D,D_std,f.diffusionFitResult.fobj.C,C_std)

%% Update lab notebook with results

% PUTING IN THE CORRECT NOTEBOOK PAGE TITLE IS VERY IMPORTANT
% --------------------
notebook = 'Matt Lab Notebook';
folder = 'Experiments';
page_title = '2025-10-01 Diffusion of CO2 in PMIM NTF2 at 80 C';
pre_image_fig_num = 1;
post_image_fig_num = 2;
comparison_image_fig_num = 5015;
uptake_curve_fig_num = 319;
diffusion_fitting_result_fig_num = 144;
water_content_fig_num = 4745;
% --------------------

obj = labarchivesCallObj('notebook',notebook,...
    'folder',folder,...
    'page',page_title);

% parameters text entry
obj.addEntry('plain text entry',...
    sprintf(['Data post-processing parameters for easy access:\n\n'...
    '* Radius used: %.3f μm\n\n'...
    '* Displacement used: %.3f μm\n\n'...
    '* Displacement defined: %s\n\n'...
    '* Measured temperature: %.2f ± %.2f ºC'],...
    f.radius, f.displacement, image_processing.pre.dx_definition,...
    f.temperature, f.temperature_std));

% pre photo
% ----------
figure(pre_image_fig_num)
% ----------
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + image_processing.pre.radius + "mm, " + "dx = " + image_processing.pre.displacement + "mm.";
obj = obj.updateFigureAttachment('caption',caption);

if post_image_complete
    caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
    caption = caption + "radius = " + image_processing.post.radius + "mm, " + "dx = " + image_processing.post.displacement + "mm.";
    
    % post photo
    % ----------
    figure(post_image_fig_num)
    % ----------
    obj = obj.updateFigureAttachment('caption',caption);
    
end

% comparison photo
% ----------
figure(comparison_image_fig_num)
% ----------
obj = obj.updateFigureAttachment('caption',"Comparison photo between the pre- and post-diffusion image. Differences are highlighted in green and pink. Images were translated and rotated to find best match using imregconfig.");

% uptake curve
% ----------
figure(uptake_curve_fig_num)
% ----------
obj = obj.updateFigureAttachment;

% single diffusion fitting result
% ----------
figure(diffusion_fitting_result_fig_num)
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
figure(water_content_fig_num)
% ----------
obj = obj.updateFigureAttachment('caption',"Water content throughout the experiment.");

%% ---- END OF STEP 3 ----

%% Step 4: Fit for double diffusion coefficient

% Proceed to the double_diffusion_fitting.m script for this part. You will
% need to have run this entire script for it to work as it will draw from
% the .mat files you saved here.