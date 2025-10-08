function concs = concOverTime(fittedSpectra,path_length,varargin)

%     path length MUST be in microns

if mod(numel(varargin),2) ~= 0
    error("Watch for unpaired name/value pairs!")
end
epsilon = 1050;
hot_band = false;
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "epsilon"
            epsilon = val;
        otherwise
            error("Invalid name/value pair")
    end
    varargin = varargin(3:end);
end
if isempty(fittedSpectra)
    error('You do not have any fitted spectra. Fit the gas lines out first.')
    return
end
n_spectra = numel(fittedSpectra);
OD = zeros(1,n_spectra);
for ii = 1:n_spectra
    temp = fittedSpectra(ii).fobj;
    fcn = co2GasLineFitFunction(fittedSpectra(ii).x,...
        temp.center,temp.w_g,temp.w_l,temp.a1,temp.a2,0,0,0);
    OD(ii) = max(fcn);
end
% molar absorptivity in units of M^-1 cm^-1
% converts path length to centimeters
concs = OD./(epsilon*path_length*1e-4);
end