%%

% input data at the start
% -------------------------------
isilon_path = "/Volumes/CHEM-SGR/"; % path to CHEM-SGR on YOUR computer
experiment_metadata.date_of_experiment = "2025-10-07";
experiment_metadata.run_number = 1; % in case of multiple runs in a day
experiment_metadata.spectra_range = 1:344;
experiment_metadata.file_prefix = "PEGDA_20251007_rt_";
experiment_metadata.temperature_log_filename = "TemperatureLog[8_24_23_AM][10_7_2025].log";
experiment_metadata.path_length = 10.2680; % μm
experiment_metadata.time_delay = 30; % seconds
experiment_metadata.sample_name = "PEGDA";
experiment_metadata.your_name = "Matt";
experiment_metadata.temperature_setpoint = NaN;
% -------------------------------

% ===================================================
%   Step 1 - Load in FTIR data
% ===================================================

% get the year for directory purposes
year_of_experiment = year(datetime(experiment_metadata.date_of_experiment));

% path to the data
data_path = isilon_path + "sgr-ftir.chem.pitt.edu/" + year_of_experiment + "/" + experiment_metadata.date_of_experiment;
try
    cd(isilon_path)
catch
    error("An error ocurred accessing the data directory. Make sure you are connected to Isilon.")
end

% load in the spectra and frequency axis
fprintf("Loading the data...\n")
cd(data_path)
tic
[FTIR_data.spectra,freq] = LoadSpectra(data_path,experiment_metadata.file_prefix,experiment_metadata.spectra_range);
T = toc;
freq = freq(:,1);
if freq(2) - freq(1) > 0
    FTIR_data.freq = flip(freq);
else
    FTIR_data.freq = freq;
end

% % initialize FTIR object
% f = FTIRexperiment(sub_data,freq,NaN,path_length,NaN,...
%     time_delay,sample_name,date_of_experiment,your_name);
% f = f.timeAxis(data_path,file_prefix,spectra_range);
% f.temperature_setpoint = temperature_setpoint;

% read the temperature
 try
    temp_log = readmatrix(experiment_metadata.temperature_log_filename);
catch
    warning("No temperature log found. A temperature was not recorded.")
    experiment_metadata.temperature = NaN;
    temperature_std = NaN;
    temp_log = [];
end
if ~isempty(temp_log)
    temp_log = temp_log(:,2);
    temp_log = temp_log(~isnan(temp_log));
    experiment_metadata.temperature = mean(temp_log);
    experiment_metadata.temperature_std = std(temp_log);
    fprintf("Temperature recorded: %.3f ± %.3f ºC\n",experiment_metadata.temperature,experiment_metadata.temperature_std)
else
    fprintf("Cannot read temperature log file until experiment is complete.\n")
end

fprintf("Time to load %i spectra: %.3f mins\n",size(FTIR_data.spectra,2),T/60);
fprintf("Successfully imported %i spectra.\n",size(FTIR_data.spectra,2))

% Results of step 1: spectra loaded and metadata organized
clearvars -except experiment_metadata FTIR_data data_path isilon_path

% ===================================================
%  Step 2 - Fit the FTIR peaks to get concentration
% ===================================================

FTIR_data.trial_spectrum = size(FTIR_data.spectra,2);
FTIR_data.fit_range = [2290 2390];

% get user input for start point parameters
start_point = true;
while start_point
    
    sp_vec = input(sprintf(['Enter start point values as a vector:\n'...
        '[center\twg\twl\ta1\ta2\ta3\tc0\tc1]\n']));
    
    sp.center = sp_vec(1);
    sp.wg = sp_vec(2);
    sp.wl = sp_vec(3);
    sp.a1 = sp_vec(4);
    sp.a2 = sp_vec(5);
    sp.a3 = sp_vec(6);
    sp.c0 = sp_vec(7);
    sp.c1 = sp_vec(8);
    
    fprintf("Chosen start point values:\n")
    sp
    
    % take the difference spectrum from spectrum 1
    s = FTIR_data.spectra(:,FTIR_data.trial_spectrum) - FTIR_data.spectra(:,1);
    
    %fit function requires fliipped inputs
    freq = flip(FTIR_data.freq);
    s = flip(s);
    
    %get x and y data for the fit
    ind1 = find(freq>=FTIR_data.fit_range(1) & freq<FTIR_data.fit_range(2));
    x = freq(ind1);
    ydata = s(ind1);
    
    %plot the fitted function using user parameters
    yfit = co2GasLineFitFunction(x,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1);
    res = ydata-yfit;
    sse = sum(res.^2);
    
    figure(182);clf
    plot(x,ydata,'o',x,yfit,x,res-0.1,'r-o')
    title(sprintf("Start point for FTIR fitting using spectrum %i",FTIR_data.trial_spectrum))
    
    % Ask the user if the start point is good
    answer = input(sprintf("Start point for FTIR fitting using spectrum %i. Retry?\t",FTIR_data.trial_spectrum),'s');
    if answer == "no"
        start_point = false;
        FTIR_data.sp = sp;
        close(figure(182))
    end
    
end

% Do the FTIR peak fit
T = tic;
diff_spectra = FTIR_data.spectra - FTIR_data.spectra(:,1);
FTIR_data.fittedSpectra = co2GasLineFit(diff_spectra,FTIR_data.freq,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1,...
    "TolFun",1e-8,"TolX",1e-8,"MaxFunEvals",5000);
stop = toc(T);

fprintf("Fitting took %.3f seconds.\n",stop)

% get time axis
time_axis = getTimeAxis(data_path,experiment_metadata.file_prefix,experiment_metadata.spectra_range);

% convert spectra to conc over time
conc_over_time = concOverTime(FTIR_data.fittedSpectra,experiment_metadata.path_length);

% Plotting the uptake curve for viewing
figure(319);clf
n_spectra = size(FTIR_data.fittedSpectra,2);
for ii = 1:n_spectra
    temp = FTIR_data.fittedSpectra(ii).fobj;
    pf = co2GasLineFitFunction(FTIR_data.fittedSpectra(ii).x,temp.center,temp.w_g,temp.w_l,...
        temp.a1,temp.a2,0,0,0);
    plot(subplot(2,1,1),FTIR_data.fittedSpectra(ii).x,pf)
    hold on
end
% fitted spectra display
title('Baseline-corrected fitted spectra with gas lines removed.')
xlabel('Wavenumbers (cm^{-1})')
ylabel('Absorbance (AU)')
box off
set(gca,'TickDir','out')
hold off
% uptake curve
plot(subplot(2,1,2),time_axis,conc_over_time,'o-','color','blue');
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

% Water content analysis

% wavenumber range of the water peaks
region_of_interest = [3150 3900];

% integration range
integration_range = [3300 3650];

% frequency point to track
freqpt = 3513;

% Plot of water over time
y_max = zeros(1,numel(FTIR_data.fittedSpectra));
y_int = y_max;
t = time_axis/3600; % convert to hours

figure(4745);clf
clear sub_freq sub_data
subplot(2,1,2)
hold on

for ii = 1:numel(t)
    
    [sub_data(:,ii),sub_freq] = getDataSubset(FTIR_data.freq,FTIR_data.spectra(:,ii),region_of_interest);
    
    sub_data(:,ii) = baselineCorrect(sub_freq,sub_data(:,ii),"value",3200);
    
    curve_color = [ii/numel(t) 0 (numel(t)-ii)/numel(t)];
    
    plot(sub_freq,sub_data(:,ii),'Color',curve_color)
    
    y_max(ii) = getDataSubset(sub_freq,sub_data(:,ii),[freqpt freqpt]);
    
    y_int(ii) = trapz(getDataSubset(sub_freq,sub_data(:,ii),integration_range));
    
end
freqline = xline(freqpt,'g--','LineWidth',2);
int_range = plot([integration_range(1) integration_range(2)],[max(sub_data(:,end)) max(sub_data(:,end))],'go--','LineWidth',2);
xlim(region_of_interest)
xlabel('Frequency (cm^{-1})')
ylabel('Difference Absorbance')
legend([freqline, int_range],{'Frequency point to track','Integration range'})
box off
hold off

subplot(2,1,1)
yyaxis left
plot(t,y_max,'bo-','MarkerFaceColor','blue','MarkerEdgeColor','blue')
ylabel('Max Absorbance (O.D.)')
yyaxis right
plot(t,y_int,'ro-','MarkerFaceColor','red','MarkerEdgeColor','red')
ylabel('Integrated Intensity (A.U.)')
xlabel('Time (hr)')
legend('Absorbance at frequency point','Integrated intensity','Location','northwest')
box off

set(gcf,'Position',[946     1   648   946])
set(gcf,'Color','white')

% Clean up the workspace
clearvars -except experiment_metadata FTIR_data time_axis conc_over_time isilon_path data_path

%%
% ===================================================
%   Step 3 - Analyze the image
% ===================================================

clear image_processing
% -------------------------
pre_image_path = "/sgr-kiralux-camera.chem.pitt.edu/2025-10-07/";
pre_image_filename = "pegda_20251007_rt_pre-diffusion.tif";
post_image_path = "";
post_image_filename = "";
image_processing.pre.dx_definition = "from edge";
image_processing.pre.radius_definition = "circle fit";
image_processing.pre.units = "um";
% -------------------------

if post_image_filename == ""
    post_image_complete = false;
else
    post_image_complete = true;
end

path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
ref_filename = "scaled_reticle01.tif";

% analyze the pre image
[image_processing.pre.radius,...
    image_processing.pre.displacement,...
    image_processing.pre.other_data] = radius_from_image(isilon_path+pre_image_path+pre_image_filename,...
    path_to_ref,ref_filename,...
    "image scale",80,"dx definition",image_processing.pre.dx_definition,...
    "radius definition",image_processing.pre.radius_definition,...
    "units",image_processing.pre.units,"flag_plot",true);

fprintf("Pre-image: r = %4f μm; dx = %4f μm\n",image_processing.pre.radius,image_processing.pre.displacement)

% package up the pre-data
image_processing.pre.path = pre_image_path + pre_image_filename;

if post_image_complete
    
    % analyze the post image
    [image_processing.post.radius,...
        image_processing.post.displacement,...
        image_processing.post.other_data] = radius_from_image(isilon_path+post_image_path+post_image_filename,...
        path_to_ref,ref_filename,...
        "image scale",40,"dx definition",image_processing.pre.dx_definition,...
        "radius definition",image_processing.pre.radius_definition,...
        "units",image_processing.pre.other_data.units,"flag_plot",true);
    
    fprintf("Post-image: r = %4f μm; dx = %4f μm\n",image_processing.post.radius,image_processing.post.displacement)
    
    % package up the post-data
    image_processing.post.path = post_image_path + post_image_filename;
    image_processing.post.units = image_processing.pre.other_data.units;
    image_processing.post.dx_definition = image_processing.pre.dx_definition;
    image_processing.post.radius_definition = image_processing.pre.radius_definition;
    
    % Show the image comparison
    pre_image = imread(isilon_path + image_processing.pre.path);
    post_image = imread(isilon_path + image_processing.post.path);
    [optimizer, metric] = imregconfig('multimodal');
    post_image_moved = imregister(post_image,pre_image,'similarity',optimizer,metric);
    figure(5015);
    imshowpair(pre_image,post_image_moved)
    
end

% set the displacement and radius
if isfield('image_processing','post')
    radius = mean([image_processing.pre.radius image_processing.post.radius]);
    displacement = mean([image_processing.pre.displacement image_processing.post.displacement]);
else
    radius = image_processing.pre.radius;
    displacement = image_processing.pre.displacement;
end

fprintf("Final values: r = %.5f μm, dx = %.5f\n",radius,displacement)

% ===================================================
%   Step 4 - Fit for diffusion coefficient
% ===================================================

diffusion_start_point = true;
while diffusion_start_point
    
    % ---- User input starting values
    dx = displacement;
    ub = [1e5 10];
    lb = [0 0];
    % ----
    
    sp_vec = input(sprintf(['Diffusion start point params as vector:\n'...
        '[D\tC]\n']));
    
    sp = sp_vec;
    
    % time axis
    t = time_axis;

    % concentration data
    y = conc_over_time;
    
    % other model parameters
    A = radius;
    nmax = 150;
    rres = 50;
    rlim = 350;
    sigma = 704;
    dy = 0;
    dx_def = image_processing.pre.dx_definition;
    
    % show the starting point plot
    figure(728);clf
    plot(t,y)
    hold on
    ymodel = diffusion_moving_beam(t,sp(1),A,sp(2),nmax,sigma,dx,dy,"rlim",rlim,'dx definition',dx_def);
    plot(t,ymodel)
    res = y(:) - ymodel(:);
    plot(t,res-0.025,'ro')
    errs(1) = sum((res).^2);
    errs(2) = std(res);
    fprintf("SSE: %.5f; Std. res: %.5f\n", errs(1), errs(2))
    
    % Ask the user if the start point is good
    answer = input(sprintf("Start point for diffusion fitting. Retry?\t"),'s');
    if answer == "no"
        diffusion_start_point = false;
        close(figure(728))
    end
    
end

% Do the uptake curve fitting
fprintf("Beginning diffusion coefficient fit.\n")
diffusionFitResult = fitDiffusionCoeff(t,y,sp,lb,ub,A,nmax,sigma,dx,rlim,dx_def);
fprintf("Diffusion coefficient fit complete.\n")

% display fit result
figure(144);clf

plot(diffusionFitResult.x,diffusionFitResult.ydata,...
    'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(diffusionFitResult.x,diffusionFitResult.yfit,...
    'red','LineWidth',1.5)
residuals = diffusionFitResult.ydata(:) - diffusionFitResult.yfit(:);
plot(diffusionFitResult.x,(residuals*1 - 0.02),'o','MarkerEdgeColor','red')
legend('Data points','Fitted curve','Location','northwest')
xlabel('Time (s)')
ylabel('Concentration (M)')
hold off

% show cfit object
diffusionFitResult.fobj;

% display values and standard deviations
ci = confint(diffusionFitResult.fobj);
D_std = (ci(2,1) - ci(1,1))/4;
C_std = (ci(2,2) - ci(1,2))/4;
fprintf("\n============\nResults\n============\n\tD = %.4f ± %.4f μm^2/s\n\tC = %.4f ± %.4f M\n\n",diffusionFitResult.fobj.D,D_std,diffusionFitResult.fobj.C,C_std)

% Clean up the workspace
clearvars -except image_processing radius displacement time_axis conc_over_time FTIR_data experiment_metadata diffusionFitResult isilon_path data_path

%% Save everything - DO NOT RUN IF ALREADY DONE

% ===================================================
%   Step 5 - Update lab notebook
% ===================================================

% PUTING IN THE CORRECT NOTEBOOK PAGE TITLE IS VERY IMPORTANT
% --------------------
notebook = 'Matt Lab Notebook';
folder = 'Experiments';
page_title = '2025-10-06 Diffusion of CO2 in PEGDA';
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
    radius, displacement, image_processing.pre.dx_definition,...
    experiment_metadata.temperature, experiment_metadata.temperature_std));

% pre photo
% ----------
figure(pre_image_fig_num)
% ----------
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + image_processing.pre.radius + "μm, " + "dx = " + image_processing.pre.displacement + "μm.";
obj = obj.updateFigureAttachment('caption',caption);

% post photo
% ----------
figure(post_image_fig_num)
% ----------
caption = "Kiralux camera photo of the sample before the diffusion annotated with calculated values: ";
caption = caption + "radius = " + image_processing.post.radius + "μm, " + "dx = " + image_processing.post.displacement + "μm.";
obj = obj.updateFigureAttachment('caption',caption);

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
coeffs = coeffnames(diffusionFitResult.fobj);
ci = confint(diffusionFitResult.fobj);
units = ["um^2/s" "M"];
if numel(units) ~= numel(coeffs)
    error("Cannot match all fitting parameters with a unit.")
end
for ii = 1:numel(coeffs)
    std_devs{ii} = (ci(2,ii) - ci(1,ii))/4;
    if isnan(std_devs{ii})
        caption = caption + coeffs{ii} + " = " + diffusionFitResult.fobj.(coeffs{ii})...
            + " " + units(ii) + ", ";
    else
        caption = caption + coeffs{ii} + " = " + diffusionFitResult.fobj.(coeffs{ii})...
            + " ± " + std_devs{ii} + " " + units(ii) + ", ";
    end
end
obj = obj.updateFigureAttachment('caption',caption);

% water content
% ----------
figure(water_content_fig_num)
% ----------
obj = obj.updateFigureAttachment('caption',"Water content throughout the experiment.");

% ===================================================
%   Step 6 - Save the data to Isilon
% ===================================================

cd(data_path)
clearvars -except image_processing radius displacement time_axis conc_over_time FTIR_data experiment_metadata diffusionFitResult isilon_path data_path
save(sprintf("%s_diffusion_analysis_%04d.mat",experiment_metadata.date_of_experiment,experiment_metadata.run_number))