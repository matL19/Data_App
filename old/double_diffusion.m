%%
% must run diffusion_fitting.m first
cd('/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2025/2025-01-29')
load 2025-01-29.mat
uh = f;
clear f
%% fit for diffusion coefficient
%get parameters ready
t = uh.timePts;
%         t = t(1:end-3);
%           t = t(1:end-15);
%         t = t-t(1);
y = uh.concOverTime;
%         y = y(4:end);
%           y = y(1:end-15);
A = uh.radius;
C = uh.finalConc;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dx = 0;
dy = 0;
sp = [102 10 0.17 0.1 0 0]; % put guess here
ub = [1e5 1e5 1e3 1e3 0.5*uh.radius 0.5*uh.radius];
lb = [0 0 0 0 0 0];

figure(728);clf
plot(t,y)
hold on
plot(t,...
    two_step_diffusion(t,sp(1),sp(2),A,sp(3),sp(4),nmax,sigma,sp(5),sp(6),dy,rlim))
plot(t,diffusion_moving_beam(t,sp(1),A,sp(3),nmax,sigma,sp(5),dy,"rlim",rlim),'red')
plot(t,diffusion_moving_beam(t,sp(2),A,sp(4),nmax,sigma,sp(6),dy,"rlim",rlim),'green')


%% Actually do the fit

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter');

ft = fittype(@(D1,D2,C1,C2,dx1,dx2,t) two_step_diffusion(t,D1,D2,A,C1,C2,nmax,sigma,dx1,dx2,dy,rlim),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D1','D2','C1','C2','dx1','dx2'},...
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

uh.diffusionFitResult = out;

%% display fit result
figure(4);clf

plot(uh.diffusionFitResult.x,uh.diffusionFitResult.ydata,...
    'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(uh.diffusionFitResult.x,uh.diffusionFitResult.yfit,...
    'red','LineWidth',1.5)
residuals = uh.diffusionFitResult.yfit - uh.diffusionFitResult.ydata(:);
plot(uh.diffusionFitResult.x,(residuals*10 - 0.02),'o','MarkerEdgeColor','red')
legend('Data points','Fitted curve','Location','northwest')
xlabel('Time (s)')
ylabel('Concentration (M)')
hold off


% get confidence intervals
ci = confint(uh.diffusionFitResult.fobj);

% readout = [string(uh.diffusionFitResult.fobj.D)]
% others = ["95% Confidence Interval is "+ci(1)+" to "+ci(2)+".";...
%     "R^2 = "+string(uh.diffusionFitResult.G.rsquare)]

uh.diffusionFitResult.fobj
%% 
save('2025-01-29_two_step_fitting.mat','uh')