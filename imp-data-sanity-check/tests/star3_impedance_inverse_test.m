%Impedance Inverse Test
%
% purpose: check that the utilities work for basic impedance problems
% 
% with the parameters below, we see:
%   - 'i' -> converges quickly, doesn't care much about parameters
%   - 'o' -> not bad. requires higher frequency than 
%             expected to get close.
%   - 'io' -> worse than 'o'. helps to allow complex lambda?

clear all; clc; close all;

path_to_ios2d = '../../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

fnamebase = '../data-out/star3_ik1_nk25_tensor_data_Impedance';
fname = [fnamebase, '.mat'];

load(fname);
warning('off')

optim_opts = [];
opts = [];
opts.verbose=true;
opts.ncoeff_impedance_mult = 0.5;
bc = [];
bc.type = 'Impedance';
bc.invtype = 'o';
optim_opts.optim_type = 'gn';
optim_opts.eps_curv = 1e-3;
optim_opts.filter_type = 'gauss-conv';
optim_opts.eps_res = 1e-7;
optim_opts.eps_upd = 1e-7;
optim_opts.maxit = 10;
opts.store_src_info = true;
opts.lambdareal=false;

if(strcmpi(bc.invtype,'i'))
    A = load(fname);
    nh = 1;
    hcoefs = zeros(2*nh+1,1);
    [src_init,varargout] = rla.update_geom(A.src_info,nh,hcoefs);
    lam = [];
    opts.ncoeff_impedance_mult = 1;    
elseif(strcmpi(bc.invtype,'o'))
    src_init = [];
else
    src_init = [];
    lam = [];
    opts.ncoeff_impedance_mult = 0.5;
    
end
    
%

[inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,lam);

fnameout = [fnamebase, '_', bc.invtype, '_invsol', char(datetime), '.mat'];                      
save(fnameout,'inv_data_all','fname','-v7.3','bc','optim_opts','opts');
rla.post_process(inv_data_all,fname);
