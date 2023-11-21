%RECONSTRUCTIONS_VARY_DELTA
%
% plot the reconstructions obtained for different impedance models
%
% test_id = 12 standard dissipation plane 2 (omegamax = 40)
% test_id = 3 standard dissipation plane 1 (omegamax = 30)
% test_id = 20 standard dissipation starfish (omegamax = 10)

path_to_ios2d = '../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

clearvars 

image_to_make = 4212;

findsigma = false;
mult_epscurv_runs = false;
findfourier = false;
findneumann = false;
findconstfirst = false;

switch image_to_make    
    case 421
        test_range = [61,62,63];
        delta_list = {'$\delta = \delta_0/16$','$\delta = \delta_0/64$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_tens';
        omega_list = [5,10,40];
    
    case 4212
        test_range = [59,60,61,62];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/4$','$\delta = \delta_0/16$','$\delta= \delta_0/64$'};
        fsavebase = 'plane2_tens';
        omega_list = [5,10,40];
    
    case 4242
        test_range = [59,61,63];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_tens_constmodel';
        omega_list = [5,10,40];
        findfourier = true;

    
    case 424
        test_range = [65,59,61];
        delta_list = {'$\delta = 16\delta_0$','$\delta = \delta_0$','$\delta = \delta_0/16$'};
        fsavebase = 'plane2_tens_neumann';
        omega_list = [5,10,40];
        findneumann = true;

    case 422
        test_range = [61,62,63];
        delta_list = {'$\delta = \delta_0/16$','$\delta = \delta_0/64$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_tens_sigma_1em1';
        omega_list = [5,10,40];
        findsigma = true;
        sigmafind = 1e-1;

    case 4222
        test_range = [59,61,63];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_tens_sigma_1em1';
        omega_list = [5,10,40];
        findsigma = true;
        sigmafind = 1e-1;
    
    case 423
        test_range = [66,68,70];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_pio8';
        omega_list = [5,20,40];
        mult_epscurv_runs = true;
        epsmake = 1e-1;

    case 4232
        test_range = [66,68,70];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_pio8';
        omega_list = [5,20,40];
        mult_epscurv_runs = true;
        epsmake = 1e-1;
        findconstfirst = true;

end


fnames_abv = {};
fnames_orig = {};

for j = 1:length(test_range)
    test_id = test_range(j);
    
    if findfourier
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*fourier*.mat',test_id);
    elseif findneumann
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*Neumann*.mat',test_id);
    elseif findconstfirst
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*constfirst*.mat',test_id);
    else
        wildcardstr = sprintf('../trans-data/data-out/test_%03d*antbar*.mat',test_id);
    end
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for abv model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st(1).name,'.mat')];
    fnametmp = [fnamebase, '.mat'];

    if mult_epscurv_runs
        for jj = 1:length(st)
            fnamebase = ['../trans-data/data-out/',erase(st(jj).name,'.mat')];
            fnametmp = [fnamebase, '.mat'];
            Atmp = load(fnametmp);
            if Atmp.optim_opts.eps_curv == epsmake
                break
            end
        end
    end


    if findsigma
        for jj = 1:length(st)
            fnamebase = ['../trans-data/data-out/',erase(st(jj).name,'.mat')];
            fnametmp = [fnamebase, '.mat'];
            Atmp = load(fnametmp);
            if isfield(Atmp,'sigma') && Atmp.sigma == sigmafind
                fprintf('found!\n')
                break
            end
        end
    end
       

    fnametmp
    fnames_abv{j} = fnametmp;

    wildcardstr = sprintf('../trans-data/data-out/test_%03d*trans.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d. not generated yet?')
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st.name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_orig{j} = fnametmp;
end

strues = {};

splots = cell(length(test_range),length(omega_list));

lamcfs_all = cell(length(test_range),1);
transparams_all = cell(length(test_range),1);
kinfo_all =  cell(length(test_range),1);

for j = 1:length(test_range)
    fnametmp = fnames_orig{j};  
    load(fnametmp,'transparams_use','kinfo_use','src_info');
    omegas = (kinfo_use.k1:kinfo_use.dk:(kinfo_use.k1+(kinfo_use.nk-1)*kinfo_use.dk))/transparams_use.c2;

    strues{j} = src_info;
    transparams_all{j} = transparams_use;
    kinfo_all{j} = kinfo_use;

    fname_inv = fnames_abv{j};
    load(fname_inv,'inv_data_all');
    inv_tmp = cell2mat(inv_data_all);
    stmp1 = vertcat(inv_tmp(:).src_info_opt);

    if ~findfourier && ~findneumann
        lamcfs_all{j} = zeros(3,length(stmp1));
        for i = 1:length(stmp1)
            lamcfs_all{j}(:,i) = stmp1(i).lamcfs;
        end
    end
    for i = 1:length(omega_list)
        [~,ii] = min(abs(omegas - omega_list(i)));
        splots{j,i} = stmp1(ii);
    end

end

%%

set(0,'defaultTextInterpreter','latex');

fac = 1.1;

fig = figure(1);
clf;
set(gcf,'Position',[0,0,1200,600])
tiledlayout(1,3,'TileSpacing','Tight');

linestyles = {'b--','b:','b-.'};


st = strues{1};
xt = st.xs; yt = st.ys;

xmax = max(xt); xmin = min(xt); xc = (xmax+xmin)/2;
dx = (xmax-xmin)/2; 
x1 = xc-dx*fac; x2 = xc+dx*fac;
ymax = max(yt); ymin = min(yt); yc = (ymax+ymin)/2;
dy = (ymax-ymin)/2; 
y1 = yc-dy*fac; y2 = yc+dy*fac;

for j = 1:length(test_range)
    t = nexttile;
    plot(xt,yt,'k-')
    hold on
    legnames = {'obstacle'};
    for i = 1:3
        sp = splots{j,i};
        xs = sp.xs; ys = sp.ys;

        plot(xs,ys,linestyles{i},'LineWidth',2)
        
        legnames{i+1} = strcat("\omega = ",sprintf("%d",omega_list(i)));
    end
    
    if (j == 1)
        legend(legnames{:})
    end

    set(gca,'XTick',[],'YTick',[]);
    h = gca;
    axis equal tight
    xlim([x1,x2]); ylim([y1,y2]);
    fontsize(gca, scale=1.5);    
    title(delta_list{j})    

end

saveas(fig,[fsavebase, '_vary_delta.epsc']);
%%

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLineLineWidth',2)

fig2 = figure(2);
clf;
set(gcf,"Position",[0,0,800,600])
t = tiledlayout(2,2,'TileSpacing','Tight');
fs = 16;

%title(t,'recovered parameters','FontSize',fs,'Interpreter','latex')

linestyles = {'b-','b:','b-.','b--'};

dmax = 0;
cmax = 0;
rmax = 0;
crmax = 0;

for j = 1:length(test_range)

    kinfo = kinfo_all{j};
    transparams = transparams_all{j};

    delta = transparams.delta;
    rhor = transparams.rho1/transparams.rho2;
    cr = transparams.c1/transparams.c2;
    
    k1 = kinfo.k1;
    dk = kinfo.dk;
    nk = kinfo.nk;

    kh = k1:dk:(k1+(nk-1)*dk);
    omegas = kh*transparams.c2;

    lamcfs = lamcfs_all{j};
    deltahats = zeros(nk,1);
    rhorhats = zeros(nk,1);
    crhats = zeros(nk,1);    

    for i = 1:length(omegas)
        [deltahats(i),rhorhats(i),crhats(i)] = convert_lamabv_beta_to_phys(lamcfs(:,i),omegas(i));
    end

    dmax = max(dmax,max(deltahats));
    rmax = max(rmax,max(rhorhats));
    cmax = max(cmax,max(crhats));
    crmax = max(crmax,max(crhats.*rhorhats));

    t = nexttile(1);

    plot(omegas,deltahats,linestyles{j})
    hold on
    plot(omegas,delta*ones(size(omegas)),'k-')
    ylim([0,dmax])
    set(gca, 'YScale', 'log')

    ylabel('$\delta$','FontSize',fs)
    
    t = nexttile(2);
    
    plot(omegas,rhorhats,linestyles{j})
    hold on
    plot(omegas,rhor*ones(size(omegas)),'k-')
    ylim([0,rmax])
    ylabel('$\rho_r$','FontSize',fs)
    
    t = nexttile(3);
    
    plot(omegas,crhats,linestyles{j})
    hold on
    plot(omegas,cr*ones(size(omegas)),'k-')
    ylim([0,cmax])
    ylabel('$c_r$','FontSize',fs)
    xlabel('$\omega$','FontSize',fs)

    t = nexttile(4);

    if j == 1
        plot(omegas,rhor*cr*ones(size(omegas)),'k-')
        hold on
    end

    plot(omegas,rhorhats.*crhats,linestyles{j})
    
    ylim([0,crmax])
    ylabel('$\rho_r c_r$','FontSize',fs)
    xlabel('$\omega$','FontSize',fs)

    if (j == length(test_range))
        
        legend('actual value',delta_list{:},'Interpreter','latex')
    end    
    

end

saveas(fig2,[fsavebase, '_imp_params.epsc']);
