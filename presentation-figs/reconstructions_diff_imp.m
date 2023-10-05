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

test_range = [20,3,12];

geo_names = {'star','plane1','plane2'};
imp_list = {'Fourier','$\lambda_{\textrm{CH}}$','$\lambda_{\textrm{ABV}}$'};

fnames_f = {};
fnames_abv = {};
fnames_ck = {};
fnames_orig = {};

for j = 1:length(test_range)
    test_id = test_range(j);
    
    wildcardstr = sprintf('../trans-data/data-out/test_%03d*antbar*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for abv model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st.name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_abv{j} = fnametmp;

    wildcardstr = sprintf('../trans-data/data-out/test_%03d*constkappa*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for ck model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st.name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_ck{j} = fnametmp;

    wildcardstr = sprintf('../trans-data/data-out/test_%03d*fourier*.mat',test_id);
    st = dir(wildcardstr);
    if isempty(st)
        warning('no output file found matching test id %d for fourier model. not generated yet?',test_id)
        return
    end
    fnamebase = ['../trans-data/data-out/',erase(st.name,'.mat')];
    fnametmp = [fnamebase, '.mat'];
    fnames_f{j} = fnametmp;

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
splots = {};

omega_lists = cell(3,1);
omega_lists{1} = [2,5,10];
omega_lists{2} = [5,10,20];
omega_lists{3} = [5,10,20];

for j = 1:length(test_range)
    fnametmp = fnames_orig{j};
    load(fnametmp,'transparams_use','kinfo_use','src_info');
    omegas = (kinfo_use.k1:kinfo_use.dk:(kinfo_use.k1+(kinfo_use.nk-1)*kinfo_use.dk))/transparams_use.c2;

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

fac = 1.3;

for j = 1:length(test_range)
    fig = figure(j);
    clf;

    tiledlayout(3,3,'TileSpacing','Compact');
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
        for l = 1:3
            sp = splots{j}{i,l};
            xs = sp.xs; ys = sp.ys;

            t = nexttile;

            plot(xt,yt,'k-')
            hold on
            plot(xs,ys,'b-')
            xlim([x1,x2]); ylim([y1,y2]);
            set(gca,'XTick',[],'YTick',[]);
            h = gca;
            axis equal tight
            if (i == 1)
                title(['$\omega =$ ', sprintf('%d',omega_list(l))]);
            end
            if (l == 1)
                ylabel(imp_list{i})
            end

            fontsize(gca, scale=1.5);
            
        end
    end

    saveas(fig,[geo_names{j}, '_diff_imp_models.pdf']);
end