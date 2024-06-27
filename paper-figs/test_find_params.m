%TEST_MULTIFREQ_PARAMS
%
%

clc;
clearvars;

A = load('../trans-data/data-out/test_059_tensdata_trans.mat');
B = load('../trans-data/data-out/test_059_tensdata_trans_io_antbar3_phaseoff_20-Sep-2023 19:55:51.mat');

u_meas = A.u_meas;
nfreq = length(u_meas);
src_info = B.inv_data_all{end}.src_info_all{end};
betas = src_info.lamcfs;
omega = u_meas{end}.kh;
[delta,rhor,cr] = convert_beta_to_phys(betas,omega);
params_init = [delta,cr,rhor];

deltat = A.transparams_use.delta;
rhort = A.transparams_use.rho1/A.transparams_use.rho2;
crt = A.transparams_use.c1/A.transparams_use.c2;

%params_init = [deltat,crt,rhort];

params = impedance_abvmodel_physical_params_mf(u_meas,src_info,params_init,(nfreq-3):nfreq);

%bc = []; bc.type = 'i'; bc.invtype = 'i';
%optim_opts = []; optim_opts.filter_type = 'step_length';
%opts = []; opts.impedance_type = 'antbarphys';
%inv_data_out = rla.rla_inverse_solver(u_meas,bc,optim_opts,opts,src_info,params_init);




function [delta,rhor,cr] = convert_beta_to_phys(betas,omega)
%
% sqrt(1-1i*delta/omega)/(rhor*cr*sqrt(1+delta^2/omega^2)) = ...
%              betas(2)*sqrt(1-1i*betas(1))
% (1-1i*delta/omega)/(rhor*(1+delta^2/omega^2)) = betas(3)*(1-1i*betas(1))
%
% delta = omega*betas(1)
% rhor = 1/(betas(3)*(1+delta^2/omega^2))
% cr = 1/(rhor*betas(2)*sqrt(1+delta^2/omega^2)) 

delta = omega*betas(1);
rhor = 1/(betas(3)*(1+betas(1)^2));
cr = 1/(rhor*betas(2)*sqrt(1+delta^2/omega^2));

end