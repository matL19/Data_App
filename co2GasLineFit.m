function fittedSpectra = co2GasLineFit(spectra,freq,center,wg,wl,a1,a2,a3,c0,c1,varargin)

tolfun = 1e-6;
tolx = 1e-6;
maxfunevals = 600;
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "TolFun"
            tolfun = val;
        case "TolX"
            tolx = val;
        case "MaxFunEvals"
            maxfunevals = val;
        otherwise
            error("Invalid name/value pair")
    end
    varargin = varargin(3:end);
end

n_spectra = size(spectra,2); % number of columns

%initial guess from inputs do this before calling function
sp = [center,wg,wl,a1,a2,a3,c0,c1];
%upper and lower bounds
lb = [2300, 0.5, 0.5,   0, -1,   0, -10, -1];
ub = [2400, 4,   4,   100,  1, inf,  10,  1];

opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfunevals);
ft = fittype(@(center,w_g,w_l,a1,a2,a3,c0,c1,w) co2GasLineFitFunction(w,center,w_g,w_l,a1,a2,a3,c0,c1),...
    'independent',{'w'},'dependent','absorbance',...
    'coefficients',{'center','w_g','w_l','a1','a2','a3','c0','c1'},...
    'options',opts);

%clear out
out(n_spectra) = struct('x',[],'ydata',[],'yfit',[],'res',[],...
    'fobj',[],'G',[],'O',[]);

% start a timer
tic

% set the fit range
range1 = [2290 2390];

freq = flip(freq);
% fit each spectrum
for ii = 1:n_spectra
   
    s = flip(spectra(:,ii));
    
    % update the fitting region (x and y)
    ind1 = find(freq>=range1(1) & freq<range1(2));
    x = freq(ind1);
    ydata = s(ind1);
    
    % do the fit
    [fobj, G, O] = fit(x,ydata,ft);
    
    
    % get fit result for plotting
    yfit = fobj(x);
    
    % pack up the data and results
    out(ii).x = x;
    out(ii).ydata = ydata;
    out(ii).yfit = yfit;
    out(ii).res = ydata - yfit;
    out(ii).fobj = fobj;
    out(ii).G = G;
    out(ii).O = O;
    
    fprintf("Fitted spectrum %i of %i.\n",ii,n_spectra)
end

% stop the timer
toc

% check results
failed_fits = 0;
for ii = 1:n_spectra
    if out(ii).O.exitflag < 1
        warning('Spectrum %i did not converge!!! Results might not be trustworthy.',ii);
        failed_fits = failed_fits + 1;
    end
end
if failed_fits == 0
    fprintf("All fits were successful.\n")
end

fittedSpectra = out;

end