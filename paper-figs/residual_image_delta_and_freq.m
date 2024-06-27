%RESIDUAL_IMAGE_DELTA_AND_FREQ
%
% make image plot of residual as it depends on dissipation and 
% frequency 
%

clearvars;

path_to_ios2d = '../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));
path_to_chunkie = '../chunkie/';
run('../chunkie/startup.m');

mult_epscurv_runs = false;
findsigma = false;
findfourier = false;
findneumann = false;
findbest = false;
findconstmodel = false;

image_to_make = 432;

switch image_to_make
    case 421
        test_range = [59:63];
        fsavebase = 'ridf_plane2_tens';
        omegaplot = [5,10,40];
        deltaplot = 61:63;
    case 4213
        test_range = [59:63];
        fsavebase = 'ridf_plane2_fourier';
        omegaplot = [5,10,40];
        deltaplot = 61:63;
        findfourier = true;
    case 422
        test_range = [59:63];
        fsavebase = 'ridf_plane2_tens_sigma1em1';
        findsigma = true;
        sigmafind = 1e-1;
        omegaplot = [5,10,40];
        deltaplot = 61:63;
    case 423
        test_range = [66:70];
        fsavebase = 'ridf_plane2_reflect_pio8';
        mult_epscurv_runs = true;
        epsmake = 1e-1;
        omegaplot = [5,10,40];
        deltaplot = 66:2:70;
    case 4242
        test_range = [59:63];
        fsavebase = 'ridf_plane2_constmodel';
        findconstmodel = true;
        omegaplot = [5,10,40];
        deltaplot = 61:63;
    case 424
        test_range = [65,64,59:61];
        fsavebase = 'ridf_plane2_neumann';
        findneumann = true;
        omegaplot = [5,10,40];
        deltaplot = [65,59,61];
    case 4232
        test_range = [66:70];
        fsavebase = 'ridf_plane2_reflect_pio8';
        mult_epscurv_runs = true;
        epsmake = 1e-1;
        omegaplot = [5,10,40];
        deltaplot = 66:2:70;
    case 43
        test_range = [71:75];
        fsavebase = 'ridf_polygon_tens';
        mult_epscurv_runs = true;
        %epsmake = 1e-1;
        omegaplot = [5,10,40];
        deltaplot = 73:75;
        findbest = true; % looks for best run based on final residual
    case 432
        test_range = [71:75];
        fsavebase = 'ridf_polygon_constmodel';
        mult_epscurv_runs = true;
        %epsmake = 1e-1;
        omegaplot = [5,10,40];
        deltaplot = 73:75;
        findbest = true; % looks for best run based on final residual
        findconstmodel = true; %gets the constant model

    otherwise
        error("image_to_make undefined")
end

fnames_inv = {};
fnames_orig = {};

for j = 1:length(test_range)
    test_id = test_range(j);
    if findfourier
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*fourier*.mat',test_id);
    elseif findconstmodel
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*fourierconstmodel*.mat',test_id);
    elseif findneumann
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*Neumann*.mat',test_id);
    else
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*antbar*.mat',test_id);
    end
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d. not generated yet?','%d')
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st(1).name,'.mat')];
    fname_inv = [fnamebase, '.mat'];
    rescomp = Inf;
    
    if mult_epscurv_runs
        for jj = 1:length(st)
            fnamebase0 = ['../trans-data/data-out/',erase(st(jj).name,'.mat')];
            fnametmp0 = [fnamebase0, '.mat'];
            Atmp0 = load(fnametmp0);
            Atmp0.inv_data_all{end}.res_opt;
            Atmp0.optim_opts.eps_curv;
            if findbest 
                restmp = Atmp0.inv_data_all{end}.res_opt;
                if restmp < rescomp
                    rescomp = restmp;
                    Atmp = Atmp0;
                    fnamebase = fnamebase0;
                    fname_inv = fnametmp0;
                end
            elseif Atmp0.optim_opts.eps_curv == epsmake
                Atmp = Atmp0;
                fnamebase = fnamebase0;
                fname_inv = fnametmp0;
                fprintf('found desired epscurv for test %d\n',test_id)
                break
            end
        end
    end

    if findsigma
        Atmp = load(fname_inv);
        for jj = 2:length(st)
            if isfield(Atmp,'sigma') && Atmp.sigma == sigmafind
                fprintf('found!\n')
                break
            else
                fnamebase = ['../trans-data/data-out/',erase(st(jj).name,'.mat')];
                fname_inv = [fnamebase, '.mat'];
            end
        end
    end

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

deltas_to_plot = [];

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

    if contains(fsavebase,'polygon')
        % true object is saved as a chunkgraph
        src_info = [];
        load(fname_orig,'cg');
        src_info.xs = cg.r(1,:);
        src_info.ys = cg.r(2,:);
    end
    strues{j} = src_info;

    if any(test_range(j) == deltaplot)
        deltas_to_plot = [deltas_to_plot,deltas(j)];
    end

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

fig = figure(); 
clf;

set(gcf,'Position',[0,0,1200,600])

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
title('Relative Residual --- $R_\omega(\hat{\Gamma},\hat{\lambda})$')
fontsize(gca, scale=2)

nexttile

h=pcolor(oo_pad,dd_pad,errs_pad);
set(h,'EdgeColor','none')
set(gca,'YScale','log')
set(gca,'ColorScale','log')
xlabel('$\omega$')
clim([c1,c2]);
colormap(brewermap([],'YlGnBu'))
colorbar
title('Error in Obstacle --- $E(\hat{\Gamma})$') 
fontsize(gca, scale=2)

hold on

for i = 1:length(deltas_to_plot)
    dd = deltas_to_plot(i);
    for j = 1:length(omegaplot)
        ww = omegaplot(j);
        plot(ww,dd,'rx','MarkerSize',15,'LineWidth',2)
    end
end

saveas(fig,[fsavebase,'_both.epsc'])


