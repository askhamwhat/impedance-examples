%Complex Impedance Inverse Test
%
% purpose: check that the utilities work for complex valued impedance
% 
% ways to run (results for parameters below):
% -> set bc.invtype = 'i'. converges quickly
% -> set bc.invtype = 'o'. not full convergence, requires higher than 
%                             expected frequency to get close
% -> set bc.invtype = 'io'. similar to 'o'
%

fnamebase = '../data-out/star3_ik1_nk25_tensor_data_complex_Impedance';
fname = [fnamebase, '.mat'];
fnameout = [fnamebase, '_invsol', char(datetime), '.mat'];

load(fname);
addpath('../');
addpath('../../');
addpath('../../inverse-obstacle-scattering2d/');
addpath(genpath_ex('../../'));
addpath(genpath_ex('../../inverse-obstacle-scattering2d/'));
warning('off')

optim_opts = [];
opts = [];
opts.verbose=true;
opts.ncoeff_impedance_mult = 0.5;
bc = [];
bc.type = 'Impedance';
bc.invtype = 'io';
optim_opts.optim_type = 'gn';
optim_opts.eps_curv = 1e-3;
optim_opts.filter_type = 'gauss-conv';
optim_opts.eps_res = 1e-7;
optim_opts.eps_upd = 1e-7;
optim_opts.maxit = 10;
opts.store_src_info = true;
opts.lambdareal = false;

%

if(strcmpi(bc.invtype,'i'))
    A = load(fname);
    nh = 1;
    hcoefs = zeros(2*nh+1,1);
    [src_init,varargout] = rla.update_geom(A.src_info,nh,hcoefs);
    laminit = [];
elseif(strcmpi(bc.invtype,'o'))
    src_init = [];
    laminit = lam;
else
    src_init = [];
    laminit = [];
end
    


[inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,laminit);

fnameout = [fnamebase, '_', bc.invtype, '_invsol', char(datetime), '.mat'];                      
save(fnameout,'inv_data_all','fname','-v7.3','bc','optim_opts','opts');
rla.post_process(inv_data_all,fname);


