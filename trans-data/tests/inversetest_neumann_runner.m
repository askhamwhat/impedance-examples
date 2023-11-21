function fnameout = inversetest_neumann_runner(test_id,invtype,eps_curv,sigma)
%INVERSETEST_NEUMANN_RUNNER a relatively stable selection of optimization
% parameters
%
% sigma - noise level (default 0)

if nargin < 2 || isempty(invtype)
    invtype = 'io';
end
if nargin < 3 || isempty(eps_curv)
    eps_curv = 1e-1;
end
if nargin < 4 || isempty(sigma)
    sigma = 0;
end

% select data set to load and some overrides...

% load everything

path_to_ios2d = '../../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

wildcardstr = sprintf('../data-out/test_%03d*trans.mat',test_id);
st = dir(wildcardstr);
if isempty(st)
    warning('no output file found matching test id %d. not generated yet?')
    return
end
fnamebase = ['../data-out/',erase(st.name,'.mat')];
fname = [fnamebase, '.mat'];

A = load(fname);

% corrupt data with noise (if requested)

u_meas = A.u_meas;
for j = 1:length(u_meas)
    utmp = u_meas{j}.uscat_tgt;
    scal = max(abs(utmp(:)));
    u_meas{j}.uscat_tgt = utmp + sigma*scal*randn(size(utmp),"like",utmp);
end

% set some reasonable options

optim_opts = [];
opts = [];
opts.verbose=true;

rho1 = A.transparams_use.rho1; rho2 = A.transparams_use.rho2;


bc = [];
bc.type = 'Neumann';
bc.invtype = invtype;
optim_opts.optim_type = 'min(gn,sd)';
optim_opts.eps_curv = eps_curv;
optim_opts.filter_type = 'min(step_length,gauss-conv)';
optim_opts.eps_res = 1e-4;
optim_opts.eps_upd = 1e-4;
optim_opts.maxit = 40;
optim_opts.stepfac = 8;
optim_opts.maxit_filter = 3;
opts.store_src_info = true;
opts.constphasefactor = false;

% modify some options for certain tests 

% % nothing yet

% GO!

if (strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io') || strcmpi(bc.invtype,'io-vp'))
    src_init = [];
elseif strcmpi(bc.invtype,'i')
    rlam = 2*pi/real(A.u_meas{end}.kh);
    opts_geo = []; opts_geo.nppw = 8;
    opts_geo.rlam = rlam;
    length(A.src_info.xs)
    src_init = rla.update_geom(A.src_info,0,0,opts_geo);
    length(src_init.xs)
else
    src_init = [];
end
    
%

[inv_data_all,~] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init);

sigstring = "sigma" + sprintf("%.1e",sigma);
fnameout = fnamebase + "_Neumann_" + bc.invtype + "_" + sigstring + "_" + string(datetime) + ".mat";
save(fnameout,'inv_data_all','fname','-v7.3','bc','optim_opts','opts','sigma');

end