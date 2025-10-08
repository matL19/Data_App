% This is a guided script meant to walk you through the process of fitting
% for a diffusion coefficient using object oriented programming
%% First, load in some already analzyed data
cd('/Volumes/CHEM-SGR/sgr-ftir.chem.pitt.edu/2024/2024-07-09')
load 2024-07-09.mat

%% fit for diffusion coefficient
%get parameters ready
% f.timeAxis;
t = f.timePts;
y = f.concOverTime;

options.D = fitParamBnd('D',55.8,0,1e4,'free');
options.gel_radius = fitParamBnd('gel_radius',f.radius,0,5000,'fixed');
options.final_conc = fitParamBnd('final_conc',0.281,0,1e3,'free');
options.dx = fitParamBnd('dx',0,0,0.5*f.radius,'free');
options.nmax = fitParamBnd('nmax',150,0,1000,'fixed');
options.rlim = fitParamBnd('rlim',350,0,100*f.radius,'fixed');
options.beam_size = fitParamBnd('beam_size',704,0,1e4,'fixed');
options.t = t;
options.t0 = fitParamBnd('t0',0,-1e4,1e4,'fixed');
options.dataMatrix = f.concOverTime;

fitting_obj = movingBeamMdl(options);


%%
figure(1);clf
plot(fitting_obj.t,fitting_obj.dataMatrix)
hold on
fitting_obj = fitting_obj.calcCurve([0.219 0 58]);
plot(fitting_obj.t,fitting_obj.simMatrix)
legend('real data','simulated curve','location','northwest')
str = "D = " + fitting_obj.paramStruct.D + " C = " + fitting_obj.paramStruct.final_conc ...
    + " dx = " + fitting_obj.paramStruct.dx;
title(str)

%%
fitting_obj = fitting_obj.fitFunction;