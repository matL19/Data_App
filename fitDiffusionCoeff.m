function out = fitDiffusionCoeff(t,y,sp,lb,ub,A,nmax,sigma,dx,rlim,dx_def)

ymodel = diffusion_moving_beam(t,sp(1),A,sp(2),nmax,sigma,dx,0,"rlim",rlim,'dx definition',dx_def);
res = y(:) - ymodel(:);
errs(1) = sum((res).^2);
errs(2) = std(res);

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter',...
    'Weights',ones(size(y))*1/errs(2),...
    'TolFun',1e-16,...
    'TolX',1e-16);

ft = fittype(@(D,C,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,0,"rlim",rlim,'dx definition',dx_def),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D','C'},...
    'options',opts);

%set up structure for storing output
out = struct('x',[],'ydata',[],'yfit',[],'res',[],...
    'fobj',[],'G',[],'O',[]);

tic

%do the fit
[fobj,G,O] = fit(t(:),y(:),ft);

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

end