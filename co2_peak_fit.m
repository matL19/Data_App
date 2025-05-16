function [ftir_peak_fitting,params] = co2_peak_fit(spectra,freq,range,sp)

n_spectra = size(spectra,2); % number of columns

if freq(2) - freq(1) < 0
    freq = flip(freq);
end

%upper and lower bounds
lb = [2300, 0.5, 0.5,   0, -1,   0, -10, -1];
ub = [2400, 4,   4,   100,  1, inf,  10,  1];

opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',cell2mat(struct2cell(sp)));
%     'Display','Iter');
ft = fittype(@(center,w_g,w_l,a1,a2,a3,c0,c1,w) co2GasLineFitFunction(w,center,w_g,w_l,a1,a2,a3,c0,c1),...
    'independent',{'w'},'dependent','absorbance',...
    'coefficients',{'center','w_g','w_l','a1','a2','a3','c0','c1'},...
    'options',opts);

%clear out
out(n_spectra) = struct('x',[],'ydata',[],'yfit',[],'res',[],...
    'fobj',[],'G',[],'O',[]);

% fit each spectrum
for ii = 1:n_spectra
    
    s = flip(spectra(:,ii));
    
    % update the fitting region (x and y)
    ind1 = find(freq>=range(1) & freq<range(2));
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
    
    if out(ii).O.exitflag < 1
        warning('Spectrum %i did not converge!!! Results might not be trustworthy.\n',ii);
    else
        fprintf("Fitted spectrum %i of %i.\n",ii,n_spectra)
    end
    
end

clear ftir_peak_fitting params
ftir_peak_fitting = out;
params.range = range;
params.sp = sp;
params.lb = lb;
params.ub = ub;
params.opts = opts;
params.fittype = ft;

end