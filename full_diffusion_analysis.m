%% FULL DATA ANALYSIS of diffusion experiment - from start to finish -

% ---- Variables that require user manual input will be surrounded ----
% like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.

% ----
date_of_experiment = "2025-08-01";
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
pre_image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-07-30/";
% ----

% ---- name of the pre-diffusion image file (use the .tif one) ----
pre_image_filename = "polyemimntf2_ovntdry_50c_07302025.tif";
% ----

%% Process the image

% Get the radius from the image
% ------------
[pre_radius_pixels,pre_displacement_pixels,pre_other] = radius_from_image(isilon_path+pre_image_path+pre_image_filename,"image scale",80,"dx definition","from center","radius definition","centroid");
% ------------

%% Get the scaling

units = "um";
path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
pre_um_per_pixel = get_distance_per_pixel(path_to_ref,"scaled_reticle01.tif",[1476 1939;824 828],2000,units,"image scale",10);

%% Calculate the radius
pre_radius_mm = pre_radius_pixels*pre_um_per_pixel/1000;
pre_displacement_mm = pre_displacement_pixels*pre_um_per_pixel/1000;
fprintf("r = %4f mm; dx = %4f mm\n",pre_radius_mm,pre_displacement_mm)

%% Process the post image

% ---- is the run complete and you have a post image?
post_image_complete = false;
% ----

if post_image_complete
    
    % ---- path to the image within CHEM-SGR (end with a slash) ----
    post_image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-08-25/";
    % ----
    
    % ---- name of the pre-diffusion image file (use the .tif one) ----
    post_image_filename = "PMIM NTF2 w FTIR windows 20250821 rt post-diffusion.tif";
    % ----
    
    % Get the radius from the image
    % ------------
    [post_radius_pixels,post_displacement_pixels,post_other] = radius_from_image(isilon_path+post_image_path+post_image_filename,"image scale",40,"dx definition","from center","radius definition","centroid");
    % ------------
    
    % Get the scaling
    
    units = "um";
    path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
    post_um_per_pixel = get_distance_per_pixel(path_to_ref,"scaled_reticle01.tif",[1476 1939;824 828],2000,units,"image scale",10);
    
    % Calculate the radius
    post_radius_mm = post_radius_pixels*post_um_per_pixel/1000;
    post_displacement_mm = post_displacement_pixels*post_um_per_pixel/1000;
    fprintf("r = %4f mm; dx = %4f mm\n",post_radius_mm,post_displacement_mm)
    
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
image_processing.pre.radius_mm = pre_radius_mm;
image_processing.pre.displacement_mm = pre_displacement_mm;
image_processing.pre.length_per_pixel = pre_um_per_pixel;
image_processing.pre.units = units;
image_processing.pre.radius_definition = "centroid";
image_processing.pre.other_data = pre_other;

if post_image_complete
    image_processing.post.path = post_image_path + post_image_filename;
    image_processing.post.radius_mm = post_radius_mm;
    image_processing.post.displacement_mm = post_displacement_mm;
    image_processing.post.length_per_pixel = post_um_per_pixel;
    image_processing.post.units = units;
    image_processing.post.radius_definition = "centroid";
    image_processing.post.other_data = post_other;
end

cd(data_path)
save(date_of_experiment + "_image_processing_data.mat",'image_processing')
%% ---- END OF STEP 1 ----

%% Step 2: Fit the FTIR peaks to obtain the uptake curve

%% Load in the spectra
cd ~
% --- the indicies of the spectra you wish to use ----
spectra_range = [1:301];
% ----

% --- the spectra file prefix ---
file_prefix = 'ML_20250730_polEmNTF2_ovnt50c';
% ----

% --- experimental parameters ---
volume = NaN;  % in microliters
spacer_size = 25;  % in microns
gel_radius = pre_radius_mm*1000;  % as previously calculated, but can be overridden
displacement = pre_displacement_mm*1000;
time_delay = 1200;  % between spectra, in seconds
sample_name = "PMIM NTF2";
your_name = "Kallol";
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
f = FTIRexperiment(sub_data,freq,volume,spacer_size,gel_radius,...
    time_delay,sample_name,date_of_experiment,your_name);
f = f.timeAxis(data_path,file_prefix,spectra_range);
f.displacement = displacement;

fprintf("Successfully imported " + size(f.data,2) + " spectra.\n")

%% Guesses for FTIR peak fitting, by eye
% ---- Which spectrum will you match to? Usually the last one is good.
trial_spectrum = 300;
% ----

% set the fit range. Usually doesn't need to be changed
range1 = [2290 2390];

clear sp
% ---- User-input starting point values ----
sp.center = 2342;
sp.wg = 1.7;
sp.wl = 1.7;
sp.a1 = 3.8;  % main peak height
sp.a2 = 0.07; % expected Boltzmann factor for bend
sp.a3 = 0.0; % gas lines
sp.c0 = 1;
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

%% Water content analysis

% ------- wavenumber range of the water peaks (shouldn't ever change but you can)
range1 = [3350 3800];
% -------

% baseline correct using a specific freq point
% (because the water line can be negative so just taking the min isn't good
% enough)
baseline_correction_point = 3350;
% ------
dw = abs(f.freqAxis(2)-f.freqAxis(1));
freq_logical = f.freqAxis < baseline_correction_point + dw/2 & f.freqAxis > baseline_correction_point - dw/2;

figure(4744);clf
hold on
for ii = 1:size(f.data,2)
    % baseline correct with the first data point
    y = f.data(:,ii);
    y_correction = y(freq_logical);
    plot(f.freqAxis,f.data(:,ii) - y_correction)
end
xlim(range1)
xlabel('Frequency (cm^{-1})')
ylabel('Absorbance')
set(gca,'FontSize',12)
box off
set(gca,'TickDir','out')
set(gcf,'Position',[754   428   840   449])
set(gcf,'Color','white')

% Plot of water over time

% ------ integration range
range1 = [3500 3750];
% ------

% ------ frequency point to track
freqpt = 3632.08;
% ------

y_max = zeros(1,numel(f.timePts));
y_int = y_max;
t = f.timePts/3600; % convert to hours
dw = abs(f.freqAxis(2)-f.freqAxis(1));
for ii = 1:numel(t)
    freq = f.freqAxis;
    %  baseline correction
    freq_logical = f.freqAxis < baseline_correction_point + dw/2 & f.freqAxis > baseline_correction_point - dw/2;
    y = f.data(:,ii);
    y_correction = y(freq_logical);
    
    freq_logical = f.freqAxis < freqpt + dw/2 & f.freqAxis > freqpt - dw/2;
    sub_data = f.data(:,ii) - y_correction;
    sub_data = sub_data(freq_logical);
    y_max(ii) = max(sub_data);
    
    freq_logical = f.freqAxis < range1(2) & f.freqAxis > range1(1);
    sub_data = f.data(:,ii) - y_correction;
    sub_data = sub_data(freq_logical);
    y_int(ii) = trapz(sub_data);
end

figure(4745);clf

subplot(2,1,1)
hold on
yyaxis left
plot(t,y_max,'-o','MarkerFaceColor','blue','Color','blue')
ylabel('max value (O.D.)')
yyaxis right
plot(t,y_int,'-o','MarkerFaceColor','red','Color','red')
ylabel('integrated water content (A.U.)')

xlabel('time (hr)')
set(gca,'FontSize',12)
legend('Max Absorbance','Integrated Intensity','Location','northeast')

subplot(2,1,2)
hold on
curve_colors = zeros(numel(f.timePts),3);
for ii = 1:size(f.data,2)
    curve_colors(ii,:) = [(ii-1)/size(f.data,2) 0 1-(ii-1)/size(f.data,2)];
    y = f.data(:,ii);
    freq_logical = f.freqAxis < baseline_correction_point + dw/2 & f.freqAxis > baseline_correction_point - dw/2;
    y_correction = y(freq_logical);
    plot(f.freqAxis,f.data(:,ii) - y_correction,'Color',curve_colors(ii,:))
end
point = xline(freqpt,'Color','green','LineStyle',':','LineWidth',3);
int_range = plot([range1(1) range1(2)],[0.05 0.05],'Color','green','LineStyle','--','Marker','o','LineWidth',1.5);
expand = 100;
xlim([range1(1)-expand range1(2)+expand])
xlabel('frequency (cm^{-1})')
ylabel('baseline-corrected absorbance')
set(gca,'FontSize',12)
box off
set(gca,'TickDir','out')
legend([point,int_range],{'max absorbance pt','integration range'},'Location','southwest')

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
t = t(1:300);
y = f.concOverTime;
y = y(1:300);
A = f.radius;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dy = 0;

% ---- User input starting values
dx = f.displacement;  % from the image analysis. can be overridden
%     D      C
sp = [3.7   0.306]; % put guess here
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
page_title = '2025-08-01 Diffusion of CO2 in PMIM NTF2';
% ----

obj = labarchivesCallObj('notebook',notebook,...
    'folder',folder,...
    'page',page_title);

% pre photo
% ----------
figure(1)
% ----------
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + pre_radius_mm + "mm, " + "dx = " + pre_displacement_mm + "mm.";
obj = obj.updateFigureAttachment('caption',caption);

if post_image_complete
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + post_radius_mm + "mm, " + "dx = " + post_displacement_mm + "mm.";
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
            + " ± " + std_devs{ii} + " " + units(ii) + ", ";
    else
        caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
            + " " + units(ii) + ", ";
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