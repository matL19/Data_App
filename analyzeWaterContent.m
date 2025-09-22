%% General script for drying analysis of any sample


%% Samples

sample_name = "PMNTF2 EMIM, 2D (ID ac75-np07-kg84)";

filenames = [
    % date folder in Isilon     % filename 
%     "2025-09-03"                "PMNTF2 EMIM 2D windows drying 20250903T0626.SPA"
%     "2025-09-04"                "PMNTF2 EMIM 2D windows_0001.SPA"
    "2025-09-05"                "PMNTF2 EMIM 2D windows drying 20250905T0926.SPA"
    "2025-09-05"                "PMNTF2 EMIM 2D windows drying 20250905T1517.SPA"
];

% load each data
clear data timepts

path = "/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/";
cd(path + "2025-08-25");

for ii = 1:size(filenames,1)
    cd("../" + filenames(ii,1))
    data(:,ii) = LoadSpectra(path + filenames(ii,1), filenames(ii,2));
    fprintf("Loaded " + filenames(ii,2) + "\n")
    file_info = dir(filenames(ii,2));
    timepts(ii) = datetime(file_info.date);
end
timepts = seconds(timepts - timepts(1));

% get freq axis
fprintf("Loading in the freq axis from the last spectrum.\n")
[~, freq] = LoadSpectra(path + filenames(ii,1), filenames(ii,2));

fprintf("Done.\n")

%% Water content analysis

% frequency viewing range
% -------
region_of_interest = [2290 4000];
% -------

% baseline correct using a specific freq point
% -------
baseline_correction_point = 3900;
% -------

% frequency point to track
% -------
freqpt = 3558.44;
% -------

% integration range
% -------
integration_range = [3250 3550];
% -------

% set the curve colors. blue (oldest) to red (most recent)
curve_colors = zeros(numel(timepts),3);
for ii = 1:size(data,2)
    curve_colors(ii,:) = [(ii-1)/size(data,2) 0 1-(ii-1)/size(data,2)];
end

% Plot of water over time

% get the sub region of the data only in the region of interest
clear sub_data sub_freq
for ii = 1:size(data,2)
    
    % crop the data to only the region of interest
    [sub_data(:,ii) sub_freq] = getDataSubset(freq,data(:,ii),region_of_interest);
    
    % baseline correct
    sub_data(:,ii) = baselineCorrect(sub_data(:,ii));
    
    % normalize
    sub_data(:,ii) = sub_data(:,ii)/max(sub_data(:,ii));
end

% set up single point and integration arrays
y_max = zeros(1,numel(timepts));
y_int = y_max;
t = timepts/3600; % convert to hours

for ii = 1:numel(t)
   
    % get single point absorbance
    y_max(ii) = getDataSubset(sub_freq,sub_data(:,ii),[freqpt freqpt]);
    
    % get integrated water content
    y_int(ii) = trapz(getDataSubset(sub_freq,sub_data(:,ii),integration_range));
    
end

figure(4745);clf

% plot the single-point and integrated absorbances
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

title("Water content analysis of " + sample_name)

% plot the spectra over time
subplot(2,1,2)
hold on
for ii = 1:size(data,2)
    plot(sub_freq,sub_data(:,ii),'Color',curve_colors(ii,:))
end
point = xline(freqpt,'Color','green','LineStyle',':','LineWidth',3);
integration_line = plot([integration_range],[0.05 0.05],'Color','green','LineStyle','--','Marker','o','LineWidth',1.5);
expand = 100;
xlim(region_of_interest)
xlabel('Frequency (cm^{-1})')
ylabel('Normalized Absorbance (to the max ROI)')
set(gca,'FontSize',12)
box off
set(gca,'TickDir','out')
legend([point,integration_line],{'max absorbance pt','integration range'},'Location','northwest')

set(gcf,'Position',[946     1   648   946])
set(gcf,'Color','white')
%% Plot the differences

range1 = [2500 4000];
for ii = 1:numel(samples)
    clear sub_data
    figure(1)
    hold on
    % subset of the data
    [sub_data,sub_freq] = getDataSubset(freq,data{1,ii},range1);
    sub_data = baselineCorrect(sub_data);
    plot(sub_freq,sub_data/max(sub_data))
%     a = annotation('textbox',[0.55 0.1 0.5 0.5],'FitBoxToText','on','string',...
%         ["Sample: " + samples(ii) "Drying method: " + drying_method(ii)]);
%     xlim(range1)
    xlabel("frequency")
    ylabel("normalized absorbance")
    
%     exportgraphics(figure(ii),"/Users/matthewliberatore/Desktop/aa_temp_sgr/fig_" + string(ii) + ".pdf")
end

%% Show the CO2 in some

range1 = [2290 2390];
dw = abs(freq(2) - freq(1));
l = freq > range1(1) - dw/2 & freq < range1(2) + dw/2;
sub_freq = freq(l);
for ii = 1:numel(samples)
    clear sub_data
    figure(ii + 4);clf
    hold on
    % subset of the data
    sub_data = data{1,ii};
    sub_data = sub_data(l) - min(sub_data);
    plot(sub_freq,sub_data)
    % subset of the data
    sub_data = data{2,ii};
    sub_data = sub_data(l) - min(sub_data);
    plot(sub_freq,sub_data)
    legend('before', 'after')
    a = annotation('textbox',[0.55 0.1 0.5 0.5],'FitBoxToText','on','string',...
        ["Sample: " + samples(ii) "Drying method: " + drying_method(ii)]);
    xlim(range1)
    xlabel("frequency")
    ylabel("absorbance")
    
    exportgraphics(figure(ii + 4),"/Users/matthewliberatore/Desktop/aa_temp_sgr/fig_" + string(ii + 4) + ".pdf")
end