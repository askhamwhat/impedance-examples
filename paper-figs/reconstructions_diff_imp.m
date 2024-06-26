%RECONSTRUCTIONS_DIFF_IMP
%
% plot the reconstructions obtained for different impedance models
%
% test_id = 17 standard dissipation plane 2 (omegamax = 40)
% test_id = 18 standard dissipation plane 1 (omegamax = 30)
% test_id = 19 standard dissipation random (omegamax = 20)
% test_id = 16 standard dissipation starfish (omegamax = 30)

clearvars
path_to_ios2d = '../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

test_range = [16,17,18,19];

geo_names = {'starfish','plane1','plane2','random'};
imp_list = {'Fourier','$\lambda_{\textrm{CH}}$','$\lambda_{\textrm{ABV}}$'};

fnames_f = {};
fnames_abv = {};
fnames_ck = {};
fnames_orig = {};

for j = 1:length(test_range)
    test_id = test_range(j);
    
    wildcardstr = sprintf('../imp-data/data-out/test_%03d*antbar*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for abv model. not generated yet?',test_id)
        return
    end

    % filter out the constfirst experiment
    notfound = true;
    for jj = 1:length(st)
        if ~contains(st(jj).name,'constfirst')
            fnamebase = ['../imp-data/data-out/',erase(st(jj).name,'.mat')];
            notfound = false;
            break
        end
    end
    if notfound
        warning('no output file found matching test id %d for abv model. not generated yet?',test_id)
        return
    end
    
    fnametmp = [fnamebase, '.mat'];
    fnames_abv{j} = fnametmp;

    wildcardstr = sprintf('../imp-data/data-out/test_%03d*constkappa*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for ck model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../imp-data/data-out/',erase(st(end).name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_ck{j} = fnametmp;

    wildcardstr = sprintf('../imp-data/data-out/test_%03d*fourier*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for fourier model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../imp-data/data-out/',erase(st(end).name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_f{j} = fnametmp;

    wildcardstr = sprintf('../imp-data/data-out/test_%03d*impck.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d. not generated yet?')
        return
    end
    fnamebase = ['../imp-data/data-out/',erase(st(1).name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_orig{j} = fnametmp;
end

strues = {};
splots = {};

omega_lists = cell(4,1);
omega_lists{1} = [5,10,20];
omega_lists{2} = [5,10,20];
omega_lists{3} = [5,10,40];
omega_lists{4} = [5,10,20];

for j = 1:length(test_range)
    fnametmp = fnames_orig{j};
    load(fnametmp,'kinfo_use','src_info');
    omegas = (kinfo_use.k1:kinfo_use.dk:(kinfo_use.k1+(kinfo_use.nk-1)*kinfo_use.dk));

    strues{j} = src_info;

    fname_inv = fnames_f{j};
    load(fname_inv,'inv_data_all');
    inv_tmp = cell2mat(inv_data_all);
    stmp1 = vertcat(inv_tmp(:).src_info_opt);

    fname_inv = fnames_ck{j};
    load(fname_inv,'inv_data_all');
    inv_tmp = cell2mat(inv_data_all);
    stmp2 = vertcat(inv_tmp(:).src_info_opt);

    fname_inv = fnames_abv{j};
    load(fname_inv,'inv_data_all');
    inv_tmp = cell2mat(inv_data_all);
    stmp3 = vertcat(inv_tmp(:).src_info_opt);

    omega_list = omega_lists{j};

    splots{j} = cell(3,length(omega_list));

    for i = 1:length(omega_list)
        [~,ii] = min(abs(omegas - omega_list(i)));
        splots{j}{1,i} = stmp1(ii);
        splots{j}{2,i} = stmp2(ii);
        splots{j}{3,i} = stmp3(ii);
    end

end

%%

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

linestyles = {'r--','m:','b-.'};

fac = 1.1;

for j = 1:length(test_range)
    fig = figure(j);
    clf;
    set(gcf,'Position',[0,0,1200,600])

    tiledlayout(1,3,'TileSpacing','Tight');
    st = strues{j};
    xt = st.xs; yt = st.ys;

    xmax = max(xt); xmin = min(xt); xc = (xmax+xmin)/2;
    dx = (xmax-xmin)/2; 
    x1 = xc-dx*fac; x2 = xc+dx*fac;
    ymax = max(yt); ymin = min(yt); yc = (ymax+ymin)/2;
    dy = (ymax-ymin)/2; 
    y1 = yc-dy*fac; y2 = yc+dy*fac;

    omega_list = omega_lists{j};

    for i = 1:3

        t = nexttile;
        plot(xt,yt,'k-','LineWidth',2)
        hold on
        
        legnames = {'obstacle'};

        for l = 1:3
            sp = splots{j}{i,l};
            xs = sp.xs; ys = sp.ys;

            plot(xs,ys,linestyles{l},'LineWidth',2)
            set(gca,'XTick',[],'YTick',[]);
            h = gca;
            axis equal tight
            legnames{l+1} = strcat("$\omega = ", sprintf("%d",omega_list(l)),"$");
            
        end

        if(i == 1)
            legend(legnames{:},'Interpreter','latex')
        end
        xlim([x1,x2]); ylim([y1,y2]);

        title(imp_list{i})

        fontsize(gca,scale=1.5);

    end


    saveas(fig,[geo_names{j}, '_diff_imp_models.epsc']);
end