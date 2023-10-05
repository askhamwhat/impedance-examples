
function [fname] = generate_transmission_tensor_data(test_id,transparams, ...
    geoinfo,kinfo,path_to_ios2d,path_to_output_folder,verbose)
%GENERATE_TRANSMISSION_TENSOR_DATA
%
% This script generates the data for a transmission scattering problem,
% using the conventions of Antoine and Barucq 2005
%
% required input:
%   test_id - up to 3 digit integer, test number (goes in filename)
%
% optional (but recommended) input:
%   transparams - struct, specify speeds (c1,c2), densities (rho1,rho2),
%             and dissipation parameter delta
%
%   transparams.c1 = float (0.5)
%   transparams.c2 = float (1.0)
%   transparams.rho1 = float (1.0)
%   transparams.rho2 = float (1.0)
%   transparams.delta = float (sqrt(3.0)*max(kh)*c2)
%
%   geoinfo - struct, specify problem geometry, sensor and incident
%                direction info
%       geoinfo.name - string ('starfish'), one of starfish, simple-plane,
%                          or complicated plane
%       geoinfo.nppw - integer (20), number of points per wavelength for
%                          boundary discretization
%       geoinfo.nrecfactor - integer (10), number of receivers at each
%                          frequency is set to floor(nrecfactor*k)
%       geoinfo.nincfactor - integer (10), number of incident directions
%                          at each frequency is set to floor(nincfactor*k)
%       geoinfo.arrangement - string ('tensor')
%                   'tensor' -> each incident direction is read at every
%                   receptor
%                   'reflect' -> each incident direction is read at
%                   receptors near the angle of arrival
%                   'transmit' -> each incident direction is read at
%                   receptors opposite the angle of arrival
%       geoinfo.angle - float (2*pi) angle around incident direction in
%                        which to receive waves
%
%       geoinfo.narm - integer, (3) if starfish domain, how many arms 
%       geoinfo.amp - float, (0.3) if starfish domain, amplitude that arms
%                           deviate from rad.
%       geoinfo.rad- float, (1.0) if starfish domain, radius of circle
%
%   kinfo - struct, specify the exterior frequencies to generate data for 
%       kinfo.k1 - float (1.0), first k
%       kinfo.dk - float (0.5), spacing in k
%       kinfo.nk - integer (20), number of frequencies 
%
%   path_to_ios2d - string ('../../inverse-obstacle-scattering2d/')
%       absolute or relative location of ios2d library on system
%   path_to_output_folder - string ('../data-out/')
%       
% output:
%   fname - string of file name 


fname = '';

if (nargin < 2 || isempty(transparams))
    transparams = [];
end
if (nargin < 3 || isempty(geoinfo))
    geoinfo = [];
end
if (nargin < 4 || isempty(kinfo))
    kinfo = [];
end
if (nargin < 5 || isempty(path_to_ios2d))
    path_to_ios2d = '../../inverse-obstacle-scattering2d/';
end
if (nargin < 6 || isempty(path_to_output_folder))
    path_to_output_folder = '../data-out/';
end
if (nargin < 7 || isempty(verbose))
    verbose = false;
end

addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));


% geometry defaults
name = 'starfish';
narm = 3;
rad = 1.0;
amp = 0.3;
nrecfactor = 10;
nincfactor = 10;
nppw = 20;
arrangement = 'tensor';
angle = 2*pi;

% defaults for simple plane
nterms = 30;

if (isfield(geoinfo,'name'))
    name = geoinfo.name;
end
if (isfield(geoinfo,'narm'))
    narm = geoinfo.narm;
end
if (isfield(geoinfo,'rad'))
    rad = geoinfo.rad;
end
if (isfield(geoinfo,'amp'))
    amp = geoinfo.amp;
end
if (isfield(geoinfo,'nrecfactor'))
    nrecfactor = geoinfo.nrecfactor;
end
if (isfield(geoinfo,'nincfactor'))
    nincfactor = geoinfo.nincfactor;
end
if (isfield(geoinfo,'nppw'))
    nppw = geoinfo.nppw;
end
if (isfield(geoinfo,'nterms'))
    nterms = geoinfo.nterms;
end
if (isfield(geoinfo,'angle'))
    angle = geoinfo.angle;
end
if (isfield(geoinfo,'arrangement'))
    arrangement = geoinfo.arrangement;
end

geoinfo_use = struct('name',name,'narm',narm,'amp',amp,'rad',rad,...
    'nrecfactor',nrecfactor,'nincfactor',nincfactor,'nppw',nppw,...
    'angle',angle,'arrangement',arrangement);

fname_geo = name;
if (strcmpi(name,'starfish'))
    fname_geo = sprintf('%s_%02d_%5.2e',name,narm,amp);
end

% kinfo defaults

k1 = 1;
dk = 0.5;
nk = 20;

if isfield(kinfo,'k1')
    k1 = kinfo.k1;
end
if isfield(kinfo,'dk')
    dk = kinfo.dk;
end
if isfield(kinfo,'nk')
    nk = kinfo.nk;
end

kinfo_use = struct('k1',k1,'dk',dk,'nk',nk);

% Set of frequencies (k_{i})
kh = k1:dk:(k1+(nk-1)*dk);

% transmission parameter defaults 

c1 = 0.5;
c2 = 1.0;
rho1 = 1.0;
rho2 = 1.0;

if isfield(transparams,'c1')
    c1 = transparams.c1;
end
if isfield(transparams,'c2')
    c2 = transparams.c2;
end
if isfield(transparams,'rho1')
    rho1 = transparams.rho1;
end
if isfield(transparams,'rho2')
    rho2 = transparams.rho2;
end

delta = sqrt(3.0)*max(kh)*c2;

if isfield(transparams,'delta')
    delta = transparams.delta;
end

transparams_use = struct('c1',c1,'c2',c2,'rho1',rho1,'rho2',rho2,'delta',delta);

fname_test = sprintf('test_%03d',test_id);

%fname = [path_to_output_folder, fname_test, '_', fname_geo,'_',fname_k, ...
%    '_tensdata_imp_ck.mat'];

if strcmpi(arrangement,'tensor')
    arr_str = '_tensdata';
elseif strcmpi(arrangement,'reflect')
    arr_str = '_reflectdata';
elseif strcmpi(arrangement,'transmit')
    arr_str = '_transmitdata';
else
    warning('unknown arrangement selected, computing nothing')
    return
end
    
fname = [path_to_output_folder, fname_test, arr_str,'_trans.mat'];

% make a copy of the geometry

if strcmpi(name,'starfish')
    starcoefs = zeros(2*narm+1,1);
    starcoefs(1) = rad;
    if narm > 0
        starcoefs(narm+1) = amp;
    end

    n  = max(300,100*narm);
    src_info = geometries.starn(starcoefs,narm,n);
elseif strcmpi(name,'smooth_plane')
    n  = max(300,30*nterms);
    src_info = geometries.smooth_plane(nterms,n);
else
    warning('geoinfo.name = %s not recognized',name);
    return
end


%
bc = [];
bc.type = 'Transmission';
bc.invtype = 'o';

%convert to antoine-barucq params
omegas = kh*c2;
rhor = rho1/rho2;
cr = c1/c2;
alphas = 1.0./(rhor*(1+1i*delta./omegas));
Ns = sqrt(1+1i*delta./omegas)/cr;
antbar_params = [];
antbar_params.omegas = omegas;
antbar_params.alphas = alphas;
antbar_params.Ns = Ns;

% transmission solver params
zks = zeros(2,nk);
as = ones(2,nk);
bs = ones(2,nk);

zks(1,:) = kh.*Ns;
zks(2,:) = kh;
bs(1,:) = alphas;

save(fname,'-v7.3','src_info','transparams_use','geoinfo_use','kinfo_use');


src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;

u_meas = cell(nk,1);
L = src_info.L;

for ik=1:nk
    n = ceil(nppw*L*abs(kh(ik))/2/pi);
    n = max(n,300);
    if (strcmpi(name,'starfish'))
       n = max(n,100*narm);
       src_info = geometries.starn(starcoefs,narm,n);
    elseif (strcmpi(name,'smooth_plane'))
       n = max(n,1000);
       src_info = geometries.smooth_plane(nterms,n);
    end
    nh = 1;
    hcoefs = zeros(2*nh+1,1);
    [src_info] = rla.update_geom(src_info,nh,hcoefs);

    % set up receivers and incident directions

    n_rec  = floor(nrecfactor*kh(ik));      
    % set target locations
    %receptors (r_{\ell})
    r_tgt = 10;
    n_tgt = n_rec;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

    % Incident directions (d_{j})
    n_dir = floor(nincfactor*kh(ik));
    t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

    if strcmpi(arrangement,'tensor')
        [t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
        t_tgt_grid = t_tgt_grid(:);
        t_dir_grid = t_dir_grid(:);
        xtgt = r_tgt*cos(t_tgt_grid);
        ytgt = r_tgt*sin(t_tgt_grid);
        tgt   = [ xtgt'; ytgt'];

        sensor_info = [];
        sensor_info.tgt = tgt;
        sensor_info.t_dir = t_dir_grid;
    elseif strcmpi(arrangement,'reflect')
        t_dir_tot = [];
        tgt_tot = [];
        for idir = 1:length(t_dir)
            tt = t_dir(idir);
            diffs = min(min(abs(t_tgt-tt),abs(t_tgt-tt-2*pi)),...
                abs(t_tgt-tt+2*pi));
            itgt = diffs < angle+1e-15;
            ntgt = nnz(itgt);
            t_dir_tmp = repmat(tt,ntgt,1);
            t_dir_tot = [t_dir_tot; t_dir_tmp];
            xtgt = r_tgt*cos(t_tgt(itgt));
            ytgt = r_tgt*sin(t_tgt(itgt));
            tgt_tmp   = [ xtgt(:).'; ytgt(:).'];

            tgt_tot = [tgt_tot, tgt_tmp];

        end
        sensor_info = [];
        sensor_info.tgt = tgt_tot;
        sensor_info.t_dir = t_dir_tot;
    elseif strcmpi(arrangement,'transmit')
        t_dir_tot = [];
        tgt_tot = [];
        for idir = 1:length(t_dir)
            tt = t_dir(idir);
            diffs = min(min(min(abs(t_tgt-tt+pi),abs(t_tgt-tt-pi)),...
                abs(t_tgt-tt+3*pi)),abs(t_tgt-tt-3*pi));
            itgt = diffs < angle+1e-15;
            ntgt = nnz(itgt);
            t_dir_tmp = repmat(tt,ntgt,1);
            t_dir_tot = [t_dir_tot; t_dir_tmp];
            xtgt = r_tgt*cos(t_tgt(itgt));
            ytgt = r_tgt*sin(t_tgt(itgt));
            tgt_tmp   = [ xtgt(:).'; ytgt(:).'];

            
            tgt_tot = [tgt_tot, tgt_tmp];

        end
        sensor_info = [];
        sensor_info.tgt = tgt_tot;
        sensor_info.t_dir = t_dir_tot;
    else
        warning('unknown arrangement selected. abort')
        return
    end
    % set up and call transmission solver
    bc.transk = zks(:,ik); bc.transa = as(:,ik); bc.transb = bs(:,ik);
    [mats,erra] = rla.get_fw_mats(kh(ik),src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh(ik),src_info,mats,sensor_info,bc,opts);
    
    u_meas0 = [];
    u_meas0.kh = kh(ik);
    u_meas0.uscat_tgt = fields.uscat_tgt;
    u_meas0.tgt = sensor_info.tgt;
    u_meas0.t_dir = sensor_info.t_dir;
    u_meas0.err_est = erra;
    u_meas{ik} = u_meas0;

    if(verbose) 
        fprintf('freq %d of %d. kh = %7.4e \t err est = %7.4e\n',ik,nk,kh(ik),erra);
    end
end

save(fname,'u_meas','-append');

