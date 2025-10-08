%% This script calculates the path length using the etaloning spectrum of your sample

%%

% etaloning spectrum filepath
% -------
filepath = "/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/2025-10-07/pegda_sample_20251007_etaloning.SPA";
% -------

[path, name, ext] = fileparts(filepath);
[data,freq] = LoadSpectra(path, name + ext);

figure(1);clf
plot(freq,data)

%% Select frequency points to use

freqpts = [2764.62 3251.57];
num_fringes = 1;

path_length_mm = 10 * num_fringes / (2*abs(freqpts(2) - freqpts(1)));

fprintf("%.4f μm\n",path_length_mm*1000)