function fnameout = inversetest_runner(test_id,ifforce_constkappa,...
    ifphaseon,ifforce_fourier,invtype,iflhp,ninner,eps_curv,sigma)
%INVERSETEST_RUNNER a relatively stable selection of optimization
% parameters
%
% sigma - noise level (default 0)

if nargin < 2 || isempty(ifforce_constkappa)
    ifforce_constkappa = false;
end
if nargin < 3 || isempty(ifphaseon)
    ifphaseon = false;
end
if nargin < 4 || isempty(ifforce_fourier)
    ifforce_fourier = false;
end
if nargin < 5 || isempty(invtype)
    invtype = 'io';
end
if nargin < 6 || isempty(iflhp)
    iflhp = false;
end
if nargin < 7 || isempty(ninner)
    ninner = 5;
end
if nargin < 8 || isempty(eps_curv)
    eps_curv = 1e-1;
end
if nargin < 9 || isempty(sigma)
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

opts.lambdareal=false;

if (rho1 == rho2)
    impedance_type = 'antbar2';
else
    impedance_type = 'antbar3';
end

opts.lambdareal=true;

if ifforce_constkappa
    impedance_type = 'constkappa';
    opts.lambdareal = false;
end

if ifforce_fourier
    impedance_type = 'fourier';
end

% FORCES A MODEST BOUND ON THE MODEL COEFFICIENTS (antbar2 should be
% deprecated anyway)
if strcmpi(impedance_type,'antbar2')
    opts.lambdaprox = @(lamcfs) sign(lamcfs).*min(abs(lamcfs),100);
elseif strcmpi(impedance_type,'antbar3') % FORCES PHYSICALITY...
    opts.lambdaprox = @(lamcfs) max(lamcfs,0);
elseif strcmpi(impedance_type,'constkappa') && iflhp
    opts.lambdaprox = @(lamcfs) min(real(lamcfs),reshape([Inf,0],size(lamcfs)))+1i*min(imag(lamcfs),0);
end

bc = [];
bc.type = 'Impedance';
bc.invtype = invtype;
opts.impedance_type = impedance_type;
opts.ncoeff_impedance_mult = 0.5;
optim_opts.optim_type = 'min(gn,sd)';
optim_opts.eps_curv = eps_curv;
optim_opts.filter_type = 'min(step_length,gauss-conv)';
optim_opts.optim_type_imp = 'sd';
optim_opts.filter_type_imp = 'step_length';
if strcmpi(impedance_type,'fourier')
    optim_opts.filter_type_imp = 'min(step_length,gauss-conv)';
end
optim_opts.eps_res = 1e-4;
optim_opts.eps_upd = 1e-4;
optim_opts.maxit = 40;
optim_opts.ninner = ninner;
optim_opts.stepfac = 8;
optim_opts.maxit_filter = 3;
opts.store_src_info = true;
opts.constphasefactor = false;
if ifphaseon
    opts.constphasefactor = true;
end

% modify some options for certain tests 

% % nothing yet

% GO!

if (strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io') || strcmpi(bc.invtype,'io-vp'))
    src_init = [];
    lam_init = [];
elseif strcmpi(bc.invtype,'i')
    rlam = 2*pi/real(A.u_meas{end}.kh);
    opts_geo = []; opts_geo.nppw = 8;
    opts_geo.rlam = rlam;
    length(A.src_info.xs)
    src_init = rla.update_geom(A.src_info,0,0,opts_geo);
    length(src_init.xs)
    lam_init = [];  
else
    src_init = [];
    lam_init = A.lamcfs;
end
    
%

[inv_data_all,~] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,lam_init);

phasestr = 'phaseoff';
if (opts.constphasefactor)
    phasestr = 'phaseon';
end
sigstring = "sigma" + sprintf("%.1e",sigma);
fnameout = fnamebase + "_" + bc.invtype + "_" + opts.impedance_type + "_" + phasestr + "_" + sigstring + "_" + string(datetime) + ".mat";
save(fnameout,'inv_data_all','fname','-v7.3','bc','optim_opts','opts','sigma');

end