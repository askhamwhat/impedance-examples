%RUN_IMPEDANCE_INVERSE_TEST
%


clearvars; clc; close all;

% select data set to load

test_id = 13;

% load everything

path_to_ios2d = '../../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

wildcardstr = sprintf('../data-out/test_%03d*impck.mat',test_id);
st = dir(wildcardstr);
if isempty(st)
    warning('no output file found matching test id %d. not generated yet?')
    return
end
fnamebase = ['../data-out/',erase(st.name,'.mat')];
fname = [fnamebase, '.mat'];

load(fname);

% set some reasonable options

optim_opts = [];
opts = [];
opts.verbose=true;
opts.ncoeff_impedance_mult = 0.5;
opts.impedance_type = 'constkappa';
if exist(impedance_type,'var')
    opts.impedance_type = impedance_type;
end
bc = [];
bc.type = 'Impedance';
bc.invtype = 'io';
optim_opts.optim_type = 'min(gn,sd)';
optim_opts.eps_curv = 1e-1;
optim_opts.filter_type = 'min(step_length,gauss-conv)';
optim_opts.eps_res = 1e-4;
optim_opts.eps_upd = 1e-5;
optim_opts.maxit = 40;
opts.store_src_info = true;
opts.lambdareal=false;

% modify options for certain tests

if (test_id >= 9 && test_id <= 11)
    opts.lambdareal = true;
end
if (test_id >= 12)
    opts.lambdareal = true;
    opts.impedance_type = 'antbar2';
end

opts.constphasefactor = false;

if (strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io') || strcmpi(bc.invtype,'io-vp'))
    src_init = [];
    lam_init = [];
elseif strcmpi(bc.invtype,'i')
    src_init = src_info;
    lam_init = [];  
else
    src_init = [];
    lam_init = lamcfs;
end
    
%

[inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,lam_init);

fnameout = [fnamebase, '_', bc.invtype, '_invsol', char(datetime), '.mat'];                      
save(fnameout,'inv_data_all','fname','-v7.3','bc','optim_opts','opts');
rla.post_process(inv_data_all,fname);
