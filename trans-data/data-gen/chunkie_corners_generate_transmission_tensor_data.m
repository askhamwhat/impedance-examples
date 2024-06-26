
function [fname] = chunkie_corners_generate_transmission_tensor_data(test_id,transparams, ...
    geoinfo,kinfo,path_to_chunkie,path_to_output_folder,verbose)
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
%       geoinfo.name - string ('polygon') only doing polygons here
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
%       geoinfo.verts - 2xnverts array of vertices of polygon,
%                       counter-clockwise order 
%
%   kinfo - struct, specify the exterior frequencies to generate data for 
%       kinfo.k1 - float (1.0), first k
%       kinfo.dk - float (0.5), spacing in k
%       kinfo.nk - integer (20), number of frequencies 
%
%   path_to_chunkie - string ('../../chunkie/')
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
if (nargin < 5 || isempty(path_to_chunkie))
    path_to_chunkie = '../../chunkie/';
end
if (nargin < 6 || isempty(path_to_output_folder))
    path_to_output_folder = '../data-out/';
end
if (nargin < 7 || isempty(verbose))
    verbose = false;
end

addpath(path_to_chunkie);
run([path_to_chunkie 'startup.m'])

% geometry defaults
name = 'polygon';
nrecfactor = 10;
nincfactor = 10;
nppw = 20;
arrangement = 'tensor';
angle = 2*pi;

verts = [2 2 1 1 0 0; 0 1 1 2 2 0];

if (isfield(geoinfo,'name'))
    name = geoinfo.name;
end
if (isfield(geoinfo,'verts'))
    verts = geoinfo.verts;
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

geoinfo_use = struct('name',name,...
    'nrecfactor',nrecfactor,'nincfactor',nincfactor,'nppw',nppw,...
    'angle',angle,'arrangement',arrangement,'verts',verts);

nverts = size(verts,2);

fname_geo = name;
if strcmpi(fname_geo,'polygon')
    fname_geo = string(fname_geo) + num2str(nverts);
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

if strcmpi(name,'polygon')
    edgends = [1:nverts; 2:nverts, 1];
    cg = chunkgraph(verts,edgends);
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

save(fname,'-v7.3','cg','transparams_use','geoinfo_use','kinfo_use');


src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;

u_meas = cell(nk,1);

for ik=1:nk

    maxchunklen = 4.0/abs(kh(ik));
    opts_ref = []; opts_ref.maxchunklen = maxchunklen;
    cg1 = cg.refine(opts_ref);

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
    b1 = bs(1,ik); b2 = bs(2,ik); q = 0.5*(1/b1+1/b2);
    ks1 = zks(1,ik); ks2 = zks(2,ik); 
    ks = zks([2,1],ik); % re-order to get out-minus-in from helm2ddiff
    repin = kernel([(1/b1)*kernel('h','d',ks1), (-1/b1)*kernel('h','s',ks1)]);
    repout = kernel([(1/b2)*kernel('h','d',ks2), (-1/b2)*kernel('h','s',ks2)]);
    cf11 = [1/(q*b2), 1/(q*b1)]; cf12 = [-1/(q*b2), -1/(q*b1)];
    cf21 = [1 1]; cf22 = [-1,-1];
    repbie = kernel([ kernel('helmdiff','d',ks,cf11), kernel('helmdiff','s',ks,cf12); ...
      kernel('helmdiff','dp',ks,cf21), kernel('helmdiff','sp',ks,cf22)  ]);

    optsmat = [];
    optsmat.nsub_or_tol = 40;
    matbie = eye(2*cg1.npt) + chunkermat(cg1,repbie,optsmat);
    matbieinv = inv(matbie);

    [tgt] = unique((sensor_info.tgt).','rows');
    tgt = tgt.';
    tgtinfo = []; tgtinfo.r = tgt;

    wts2 = cg1.wts; wts2 = repmat(wts2(:).',2,1);

    mat_sol_to_receptor = repout.eval(cg1,tgtinfo)*diag(wts2(:));

    % test error 

    src_in = [];
    src_in.r = [0.5;0.5];
    src_out = [];
    src_out.r = [1.5;1.5];

    tgtin = []; tgtin.r = [0.2 0.3; 0.3 0.2];

    uin_a = chnk.helm2d.kern(ks2,src_in,cg1,'s'); % charge inside generates field outside
    dudnin_a = chnk.helm2d.kern(ks2,src_in,cg1,'sprime');
    uout_a = chnk.helm2d.kern(ks1,src_out,cg1,'s'); % charge outside generates field inside
    dudnout_a = chnk.helm2d.kern(ks1,src_out,cg1,'sprime');

    rhs_a = zeros(2*cg1.npt,1);
    rhs_a(1:2:end) = (uin_a-uout_a)/q; % out minus in, oddly enough
    rhs_a(2:2:end) = (b2*dudnin_a-b1*dudnout_a);

    sol_a = matbieinv * rhs_a;
    
    utest = mat_sol_to_receptor*sol_a;
    utest_in = repin.eval(cg1,tgtin)*(wts2(:).*sol_a);

    uex = chnk.helm2d.kern(ks2,src_in,tgtinfo,'s'); 
    erra = norm(utest(:)-uex(:))/max(norm(uex(:)),1);

    uex_in = chnk.helm2d.kern(ks1,src_out,tgtin,'s'); 
    erra = max(erra,norm(utest_in(:)-uex_in(:))/max(norm(uex_in(:)),1));
    
    mat_bd_data_to_receptor = mat_sol_to_receptor*matbieinv;
    uscat_tgt = compute_fields_loc(ks2,bs(:,ik),cg1,mat_bd_data_to_receptor,sensor_info);
    
     
    u_meas0 = [];
    u_meas0.kh = kh(ik);
    u_meas0.uscat_tgt = uscat_tgt;
    u_meas0.tgt = sensor_info.tgt;
    u_meas0.t_dir = sensor_info.t_dir;
    u_meas0.err_est = erra;
    u_meas{ik} = u_meas0;

    if(verbose) 
        fprintf('freq %d of %d. kh = %7.4e \t err est = %7.4e\n',ik,nk,kh(ik),erra);
    end
end

save(fname,'u_meas','-append');

end


function uscat_tgt = compute_fields_loc(kh,bs,cg,mat_bd_data_to_receptor,sensor_info)
%
% This subroutine computes the scattered field and it's normal derivative
% on the boundary of the obstacle and the scattered field
% at a collection of target locations due to the incident directions
% prescribed by the user, and stores both the incident and the scattered
% fields along with their normal derivatives on the boundary, and also
% the scattered field at the sensor locations
% 
% Input:

%   sensor_info - sensor information struct
%      sensor_info.tgt(2,nmeas) - xy cooordinates of sensors
%         sensor_info.tgt(1:2,i) = xy coordinates corresponding the ith
%            measurement
%      sensor_info.t_dir(nmeas) - incident directions
%         sensor_info.t_dir(i) - is the incident direction corresponding to
%            the ith measurement
%
%
      
   [t_dir_uni,~,idir] = unique(sensor_info.t_dir);
   [tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
   nt_uni = length(tgt_uni(:,1));
   induse = itgt + (idir-1)*nt_uni;
   
   
   t_dir_uni = t_dir_uni(:);
   x_dir = cos(t_dir_uni)';
   y_dir = sin(t_dir_uni)';
   n_dir = length(x_dir);

   xs = cg.r(1,:);
   ys = cg.r(2,:);
   dxs = cg.d(1,:);
   dys = cg.d(2,:);
   ds = sqrt(sum(cg.d(:,:).^2,1));
   
   fields = [];
   fields.uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
   fields.dudninc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
   
   q = 0.5*(1/bs(1) + 1/bs(2));
   n1 = numel(xs);
   bd_data = zeros(2*n1,n_dir,'like',1.0+1i);
   bd_data(1:2:end,:) = -fields.uinc/q;
   bd_data(2:2:end,:) = -bs(2)*fields.dudninc;          
   
    uscat_tgt_all = mat_bd_data_to_receptor*bd_data;
    uscat_tgt_all = uscat_tgt_all(:);
    uscat_tgt = uscat_tgt_all(induse);

   
end

