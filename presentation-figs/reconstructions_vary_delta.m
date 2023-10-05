%RECONSTRUCTIONS_DIFF_IMP
%
% plot the reconstructions obtained for different impedance models
%
% test_id = 12 standard dissipation plane 2 (omegamax = 40)
% test_id = 3 standard dissipation plane 1 (omegamax = 30)
% test_id = 20 standard dissipation starfish (omegamax = 10)

path_to_ios2d = '../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

image_to_make = 4;

switch image_to_make    
    case 1
        test_range = [12,16,19];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/128$'};
        fsavebase = 'plane2_tens';
        omega_list = [5,10,20];

    case 2
        test_range = [51,55,58];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/128$'};
        fsavebase = 'plane2_reflect_pio8';
        omega_list = [5,20,40];

    case 3
        test_range = [38,40,50];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/128$'};
        fsavebase = 'plane2_reflect_pio4';
        omega_list = [5,20,40];

    case 4
        test_range = [59,61,63];
        delta_list = {'$\delta = \delta_0$','$\delta = \delta_0/16$','$\delta = \delta_0/256$'};
        fsavebase = 'plane2_tens';
        omega_list = [5,20,40];

end


fnames_abv = {};
fnames_orig = {};

for j = 1:length(test_range)
    test_id = test_range(j);
    
    wildcardstr = sprintf('../trans-data/data-out/test_%03d*antbar*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for abv model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st(2).name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
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

    lamcfs_all{j} = zeros(3,length(stmp1));
    for i = 1:length(stmp1)
        lamcfs_all{j}(:,i) = stmp1(i).lamcfs;
    end
    
    for i = 1:length(omega_list)
        [~,ii] = min(abs(omegas - omega_list(i)));
        splots{j,i} = stmp1(ii);
    end

end

%%

set(0,'defaultTextInterpreter','latex');

fac = 1.3;

fig = figure(1);
clf;
tiledlayout(3,3,'TileSpacing','Compact');

st = strues{1};
xt = st.xs; yt = st.ys;

xmax = max(xt); xmin = min(xt); xc = (xmax+xmin)/2;
dx = (xmax-xmin)/2; 
x1 = xc-dx*fac; x2 = xc+dx*fac;
ymax = max(yt); ymin = min(yt); yc = (ymax+ymin)/2;
dy = (ymax-ymin)/2; 
y1 = yc-dy*fac; y2 = yc+dy*fac;

for j = 1:length(test_range)

    for i = 1:3
        sp = splots{j,i};
        xs = sp.xs; ys = sp.ys;

        t = nexttile;

        plot(xt,yt,'k-')
        hold on
        plot(xs,ys,'b-')
        xlim([x1,x2]); ylim([y1,y2]);
        set(gca,'XTick',[],'YTick',[]);
        h = gca;
        axis equal tight
        if (j == 1)
            title(['$\omega =$ ', sprintf('%d',omega_list(i))]);
        end
        if (i == 1)
            ylabel(delta_list{j})
        end

        fontsize(gca, scale=1.5);
        
    end
end

saveas(fig,[fsavebase, '_vary_delta.pdf']);
%%

set(0,'defaultTextInterpreter','latex');

fig2 = figure(2);
clf;
t = tiledlayout(3,3,'TileSpacing','Compact');
fs = 14;

title(t,'recovered parameters','FontSize',fs)

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
        [deltahats(i),rhorhats(i),crhats(i)] = convert_beta_to_phys(lamcfs(:,i),omegas(i));
    end

    t = nexttile;

    plot(omegas,deltahats,'k-o')
    hold on
    plot(omegas,delta*ones(size(omegas)),'k--')

    ylabel(delta_list{j},'FontSize',fs)

    if (j == 1)
        title('$\delta$','FontSize',fs)
    end
    
    t = nexttile;
    
    plot(omegas,rhorhats,'k-o')
    hold on
    plot(omegas,rhor*ones(size(omegas)),'k--')
    ylim([0,2])

    if (j == 1)
        title('$\rho_r$','FontSize',fs)
    end    
    
    t = nexttile;
    
    plot(omegas,crhats,'k-o')
    hold on
    plot(omegas,cr*ones(size(omegas)),'k--')
    ylim([0,2])

    if (j == 1)
        title('$c_r$','FontSize',fs)
    end    
    
end

saveas(fig2,[fsavebase, '_imp_params.pdf']);

function [delta,rhor,cr] = convert_beta_to_phys(betas,omega)
%
% sqrt(1-1i*delta/omega)/(rhor*cr*sqrt(1+delta^2/omega^2)) = ...
%              betas(2)*sqrt(1-1i*betas(1))
% (1-1i*delta/omega)/(rhor*(1+delta^2/omega^2)) = betas(3)*(1-1i*betas(1))
%
% delta = omega*betas(1)
% rhor = 1/(betas(3)*(1+delta^2/omega^2))
% cr = 1/(rhor*betas(2)*sqrt(1+delta^2/omega^2)) 

delta = omega*betas(1);
rhor = 1/(betas(3)*(1+betas(1)^2));
cr = 1/(rhor*betas(2)*sqrt(1+delta^2/omega^2));

end
