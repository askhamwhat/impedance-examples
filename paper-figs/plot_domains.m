%PLOT_DOMAINS
%

clearvars;

path_to_ios2d = '../inverse-obstacle-scattering2d';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

nterms = 20;
n = 1000;

src_info_plane1 = geometries.smooth_plane(nterms,n);

nterms = 70;
n = 1000;

src_info_plane2 = geometries.smooth_plane(nterms,n);

narm = 5;
amp = 0.3;
n = 1000;

nc = narm;
coefs = zeros(2*nc+1,1);
coefs(1) = 1;
coefs(nc+1) = amp;

src_info_star = geometries.starn(coefs,nc,n);


fig = figure(1);
clf
xs = src_info_plane1.xs; xs = [xs,xs(1)];
ys = src_info_plane1.ys; ys = [ys,ys(1)];
plot(xs,ys,'k-')
axis equal tight
set(gca,"Visible","off")
saveas(fig,'plane1.pdf')

fig = figure(1);
clf
xs = src_info_plane2.xs; xs = [xs,xs(1)];
ys = src_info_plane2.ys; ys = [ys,ys(1)];
plot(xs,ys,'k-')
axis equal tight
set(gca,"Visible","off")
saveas(fig,'plane2.pdf')

fig = figure(1);
clf
xs = src_info_star.xs; xs = [xs,xs(1)];
ys = src_info_star.ys; ys = [ys,ys(1)];
plot(xs,ys,'k-')
axis equal tight
set(gca,"Visible","off")
saveas(fig,'star.pdf')