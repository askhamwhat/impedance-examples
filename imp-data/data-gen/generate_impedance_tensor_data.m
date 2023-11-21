
function [fname] = generate_impedance_tensor_data(test_id,lamcfs, ...
    impedance_type,geoinfo,kinfo,path_to_ios2d,path_to_output_folder,...
    verbose)
%GENERATE_IMPEDANCE_TENSOR_DATA
%
% This script generates the data for an impedance scattering problem with a
% constant + curvature model for the impedance function
%
% required input:
%   test_id - up to 3 digit integer, test number (goes in filename)
%   lamcfs - length 1 or 2 array of numbers or cell array of function
%                     handles
%
%               lambda = cfs(1) + cfs(2)*src_info.H
%
%             or 
% 
%               lambda = cfs{1}(k) + cfs{2}(k)*src_info.H
%
% optional (but recommended) input:
%   impedance_type - string ('constkappa')
%            other options are 'antbar2' and 'antbar3'
%   geoinfo - struct, specify problem geometry, sensor and incident
%                direction info
%       geoinfo.name - string ('starfish'), one of starfish, simple-plane,
%                          or complicated plane, or random
%       geoinfo.nppw - integer (20), number of points per wavelength for
%                          boundary discretization
%       geoinfo.nrecfactor - integer (10), number of receivers at each
%                          frequency is set to floor(nrecfactor*k)
%       geoinfo.nincfactor - integer (10), number of incident directions
%                          at each frequency is set to floor(nincfactor*k)
%       geoinfo.narm - integer, (3) if starfish domain, how many arms 
%       geoinfo.nmode - integer, (5) if random domain, how many modes
%       geoinfo.amp - float, (0.3) if starfish domain, amplitude that arms
%                           deviate from rad.
%       geoinfo.rad- float, (1.0) if starfish domain, radius of circle
%       geoinfo.receptor_shape ('circle') shape where receptors are placed,
%               defaults to circle ofradius 10. 'ellipse' places them on a
%               ellipse with semi-major and semi-minor axes 16 and 8.
%
%   kinfo - struct, specify the frequencies to generate data for 
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

if (nargin < 3 || isempty(impedance_type))
    impedance_type = 'constkappa';
end
if (nargin < 4 || isempty(geoinfo))
    geoinfo = [];
end
if (nargin < 5 || isempty(kinfo))
    kinfo = [];
end
if (nargin < 6 || isempty(path_to_ios2d))
    path_to_ios2d = '../../inverse-obstacle-scattering2d/';
end
if (nargin < 7 || isempty(path_to_output_folder))
    path_to_output_folder = '../data-out/';
end
if (nargin < 8 || isempty(verbose))
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
receptor_shape = 'circle';

% defaults for simple plane
nterms = 30;

%defaults for random
nmode = 5;

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
if (isfield(geoinfo,'nmode'))
    nmode = geoinfo.nmode;
end
if (isfield(geoinfo,'receptor_shape'))
    receptor_shape = geoinfo.receptor_shape;
end

geoinfo_use = struct('name',name,'narm',narm,'amp',amp,'rad',rad,...
    'nrecfactor',nrecfactor,'nincfactor',nincfactor,'nppw',nppw,...
    'nmode',nmode,'receptor_shape',receptor_shape);

fname_geo = name;
if (strcmpi(name,'starfish'))
    fname_geo = sprintf('%s_%02d_%5.2e',name,narm,amp);
end
if (strcmpi(name,'random'))
    fname_geo = sprintf('%s_%02d_mode',name,nmode);
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

fname_k = sprintf('k1_%5.2e_dk_%5.2e_nk_%02d',k1,dk,nk);

fname_test = sprintf('test_%03d',test_id);

%fname = [path_to_output_folder, fname_test, '_', fname_geo,'_',fname_k, ...
%    '_tensdata_imp_ck.mat'];
fname = [path_to_output_folder, fname_test,'_tensdata_impck.mat'];

% make a copy of the geometry

if strcmpi(name,'starfish')
    starcoefs = zeros(2*narm+1,1);
    starcoefs(1) = rad;
    if narm > 0
        starcoefs(narm+1) = amp;
    end

    n  = max(300,100*narm);
    src_info = geometries.starn(starcoefs,narm,n);
elseif strcmpi(name,'random')
    randcoefs = 0.1*randn(2*nmode+1,1);
    randcoefs(1) = 1;

    n  = max(300,100*narm);
    src_info = geometries.starn(randcoefs,nmode,n);
elseif strcmpi(name,'smooth_plane')
    n  = max(300,30*nterms);
    src_info = geometries.smooth_plane(nterms,n);
else
    warning('geoinfo.name = %s not recognized',name);
    return
end

% Set of frequencies (k_{i})
kh = k1:dk:(k1+(nk-1)*dk);

%
bc = [];
bc.type = 'Impedance';
bc.invtype = 'o';

plot(src_info.xs,src_info.ys)
pause(1)

save(fname,'src_info','lamcfs','geoinfo_use','kinfo_use','impedance_type');


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
   elseif (strcmpi(name,'random'))
       n = max(n,100*nmode);
       src_info = geometries.starn(randcoefs,nmode,n);
   elseif (strcmpi(name,'smooth_plane'))
       n = max(n,1000);
       src_info = geometries.smooth_plane(nterms,n);
   end
   nh = 1;
   hcoefs = zeros(2*nh+1,1);
   [src_info] = rla.update_geom(src_info,nh,hcoefs);

   if (strcmpi(impedance_type,'constkappa') || ...
           strcmpi(impedance_type,'antbar2'))
       lamcfs_ik = zeros(2,1);
   elseif (strcmpi(impedance_type,'antbar3'))
       lamcfs_ik = zeros(3,1);       
   else
       error('unknown impedance_type');
   end
   
   if (isa(lamcfs,'cell'))
       for iii = 1:length(lamcfs)
           lamcfs_ik(iii) = lamcfs{iii}(kh(ik));
       end
   else
       lamcfs_ik(1:length(lamcfs)) = lamcfs;
   end

   src_info.lamcfs = lamcfs_ik;
   ckcfs = constkappa_models_convert(src_info.lamcfs,impedance_type,kh(ik));
   src_info.lambda = ckcfs(1) + ckcfs(2)*src_info.H(:);

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

    [t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
    t_tgt_grid = t_tgt_grid(:);
    t_dir_grid = t_dir_grid(:);

    if strcmpi(receptor_shape,'circle')
        xtgt = r_tgt*cos(t_tgt_grid);
        ytgt = r_tgt*sin(t_tgt_grid);
    elseif strcmpi(receptor_shape,'ellipse')
        xtgt = 8*cos(t_tgt_grid);
        ytgt = 16*sin(t_tgt_grid);
    else
        error('unknown receptor_shape %s',receptor_shape)
    end

    tgt   = [ xtgt'; ytgt'];

    sensor_info = [];
    sensor_info.tgt = tgt;
    sensor_info.t_dir = t_dir_grid;
   
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

