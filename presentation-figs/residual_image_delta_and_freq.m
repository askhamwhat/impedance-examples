%RESIDUAL_IMAGE_DELTA_AND_FREQ
%
% make image plot of residual as it depends on dissipation and 
% frequency 
%

path_to_ios2d = '../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

image_to_make = 3;
switch image_to_make
    case 1
        test_range = 12:19;
        fsavebase = 'ridf_plane2_tens';
    case 2
        test_range = [51:58];
        fsavebase = 'ridf_plane2_reflect_pio8';
    case 3
        test_range = [59:63];
        fsavebase = 'ridf_plane2_tens';
    otherwise
        test_range = [];
end

fnames_inv = {};
fnames_orig = {};

for j = 1:length(test_range)
    test_id = test_range(j);
    wildcardstr = sprintf('../trans-data/data-out/test_%03d*antbar*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d. not generated yet?','%d')
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st(1).name,'.mat')];
    fname_inv = [fnamebase, '.mat'];
    fnames_inv{j} = fname_inv;

    wildcardstr = sprintf('../trans-data/data-out/test_%03d*trans.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d. not generated yet?',test_id)
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st.name,'.mat')];
    fname_orig = [fnamebase, '.mat'];
    fnames_orig{j} = fname_orig;
end

res_opts = {};
deltas = zeros(length(test_range),1);
stmps = {};
strues = {};

for j = 1:length(test_range)
    fname_inv = fnames_inv{j};
    load(fname_inv,'inv_data_all','fname');
    inv_tmp = cell2mat(inv_data_all);
    res_opt = vertcat(inv_tmp.res_opt);
    res_opts{j} = res_opt;

    stmp = vertcat(inv_tmp(:).src_info_opt);
    stmps{j} = stmp;

    fname_orig = fnames_orig{j};
    load(fname_orig,'transparams_use','kinfo_use','src_info');
    deltas(j) = transparams_use.delta;
    omegas = (kinfo_use.k1:kinfo_use.dk:(kinfo_use.k1+(kinfo_use.nk-1)*kinfo_use.dk))/transparams_use.c2;

    strues{j} = src_info;
end

%

% compute errs 

errs = zeros(length(test_range),length(omegas));

x1 = linspace(-1.2,1.2,300);
[xx,yy] = meshgrid(x1,x1);
h = x1(2)-x1(1);

p1 = polyshape(strues{1}.xs(:),strues{1}.ys(:));

errnrm = area(p1);

%in2 = inpolygon(xx(:),yy(:),strues{1}.xs,strues{1}.ys);

%errnrm = h*h*nnz(in2);

for j = 1:length(test_range)
    fprintf('getting errs for test %d\n',test_range(j))
    stmp = stmps{j};
    for i = 1:length(stmp)
        p2 = polyshape(stmp(i).xs,stmp(i).ys);
        p3 = subtract(p1,p2);
        p4 = subtract(p2,p1); 
        aa = area(p3) + area(p4);
        
        %tic; in1 = inpolygon(xx(:),yy(:),stmp(i).xs,stmp(i).ys);
        %diff1 = and(in1,~in2);
        %diff2 = and(in2,~in1); toc
        %aa/errnrm - h*h*(nnz(diff1)+nnz(diff2))/errnrm
        errs(j,i) = aa/errnrm;
    end
end
        
%%

set(0,'defaultTextInterpreter','latex');

nfreq = length(res_opts{1});
resmat = zeros(length(test_range),nfreq);
[dd,oo] = meshgrid(omegas,deltas);

tmp = log2(deltas);
dx = tmp(2)-tmp(1);
tmp2 = zeros(length(deltas)+1,1);
tmp2(1:end-1) = tmp-dx/2;
tmp2(end) = tmp(end)+dx/2;
deltas_pad = 2.^tmp2;

omegas_pad = zeros(length(omegas)+1,1);
dx = omegas(2)-omegas(1);
omegas_pad(1:end-1) = omegas-dx/2;
omegas_pad(end) = omegas(end)+dx/2;

[oo_pad,dd_pad] = meshgrid(omegas_pad,deltas_pad);


for j = 1:length(test_range)
    resmat(j,:) = res_opts{j};
end

resmat_pad = zeros(size(resmat)+[1 1]);
resmat_pad(1:end-1,1:end-1) = resmat;

errs_pad = zeros(size(errs)+[1 1]);
errs_pad(1:end-1,1:end-1) = errs;

c1 = min(min(errs(:)),min(resmat(:))); c1 = 5*1e-3;
c2 = max(max(errs(:)),max(resmat(:))); c2 = 1;

fig = figure(1); 
clf;

tiledlayout(1,2,"TileSpacing","tight")

nexttile 

h=pcolor(oo_pad,dd_pad,resmat_pad);
set(h,'EdgeColor','none')
set(gca,'YScale','log')
set(gca,'ColorScale','log')
xlabel('$\omega$')
ylabel('$\delta$')
clim([c1,c2]);
colormap(brewermap([],'YlGnBu'))
title('Relative Residual of Fit')
fontsize(gca, scale=1.5)

nexttile

h=pcolor(oo_pad,dd_pad,errs_pad);
set(h,'EdgeColor','none')
set(gca,'YScale','log')
set(gca,'ColorScale','log')
xlabel('$\omega$')
clim([c1,c2]);
colormap(brewermap([],'YlGnBu'))
colorbar
title('Error in Recovered Obstacle')
fontsize(gca, scale=1.5)


saveas(fig,[fsavebase,'_both.pdf'])


