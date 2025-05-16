function concs = conc_over_time(data_struct,varargin)
if mod(numel(varargin),2) ~= 0
    error("Watch for unpaired name/value pairs!")
end
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "epsilon"
            epsilon = val;
        case "hot band"
            hot_band = val;
        otherwise
            error("Invalid name/value pair")
    end
    varargin = varargin(3:end);
end
if ~exist('epsilon','var')
    epsilon = 1050;
end
if ~exist('hot_band','var')
    hot_band = false;
end
if ~isfield(data_struct,'ftir_peak_fit')
    error('You do not have any fitted spectra. Fit the gas lines out first.')
end
n_spectra = numel(data_struct.ftir_peak_fit);
OD = zeros(1,n_spectra);
for ii = 1:n_spectra
    temp = data_struct.ftir_peak_fit(ii).fobj;
    x = data_struct.ftir_peak_fit(ii).x;
    fcn = co2GasLineFitFunction(x,...
        temp.center,temp.w_g,temp.w_l,temp.a1,temp.a2,0,0,0);
    if isnumeric(hot_band)
        hb_region = [2325 2333];
        new_fcn = fcn(x > hb_region(1) & x < hb_region(2));
        OD(ii) = max(new_fcn)/hot_band;
    elseif (isstring(hot_band) || ischar(hot_band)) && hot_band == "manual"
        if max(fcn) > 0.9
            hb_region = [2325 2333];
            new_fcn = fcn(x > hb_region(1) & x < hb_region(2));
            OD(ii) = max(new_fcn)/0.07;
        else
            OD(ii) = max(fcn);
        end
    else
        OD(ii) = max(fcn);
    end
end
% molar absorptivity in units of M^-1 cm^-1
% converts path length to centimeters
concs = OD./(epsilon*data_struct.path_length*1e-4);
end