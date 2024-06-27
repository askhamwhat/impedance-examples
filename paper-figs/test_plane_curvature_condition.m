
% test if plane satisfies curvature condition
%

path_to_ios2d = '../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

npts = 1000;
nterms = 70;
src_info = geometries.smooth_plane(nterms,npts);

opts = [];
opts.n_curv = 215;
opts.eps_curv = 0.1;

[src_out,ier] = rla.update_geom(src_info,0,0,opts);

ier