%SCATTERING_PLOTS_TRANS
%
% this script generates plots of scattering problems for a penetrable
% medium with dissipation, using the full transmission model 
%
%
%

clearvars;

path_to_ios2d = '../inverse-obstacle-scattering2d';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));


% geometry settings
nppw = 20;
nterms = 70;

% frequency

khs = [30];
khmax = 40;

% transmission parameter values 

c1 = 0.5;
c2 = 1.0;
rho1 = 0.7;
rho2 = 1.2;

delta_div = 1;
delta = sqrt(3.0)*khmax*c2/delta_div;

fnamebase = sprintf('trans_plot_nterms_%d_deltadiv_%d_kh_',nterms,delta_div);

%
bctrans = [];
bctrans.type = 'Transmission';
bctrans.invtype = 'o';
bcimp = [];
bcimp.type = 'Impedance';
bcimp.invtype = 'o';

%convert to antoine-barucq params
omegas = khs*c2;
rhor = rho1/rho2;
cr = c1/c2;
alphas = 1.0./(rhor*(1+1i*delta./omegas));
Ns = sqrt(1+1i*delta./omegas)/cr;
antbar_params = [];
antbar_params.omegas = omegas;
antbar_params.alphas = alphas;
antbar_params.Ns = Ns;

% transmission solver params
nk = length(khs);
zks = zeros(2,nk);
as = ones(2,nk);
bs = ones(2,nk);

zks(1,:) = khs.*Ns;
zks(2,:) = khs;
bs(1,:) = alphas;

src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;

u_meas = cell(nk,1);

src_info = geometries.smooth_plane(nterms,1000); L = src_info.L;

x1 = linspace(-3,3,200);
[xx,yy] = meshgrid(x1,x1);
tgt = [xx(:).'; yy(:).'];

for ik=1:nk
    kh = khs(ik);
    n = ceil(nppw*L*abs(kh)/2/pi);
    n = max(n,300);
    %n = max(n,1000);
    src_info = geometries.smooth_plane(nterms,n);

    in = inpolygon(xx(:),yy(:),src_info.xs(:),src_info.ys(:));

    nh = 1;
    hcoefs = zeros(2*nh+1,1);
    [src_info] = rla.update_geom(src_info,nh,hcoefs);

    % set up receivers and incident directions
    sensor_info = [];
    sensor_info.tgt = zeros(2,0);
   
    % set up transmission matrices
    bctrans.transk = zks(:,ik); bctrans.transa = as(:,ik); bctrans.transb = bs(:,ik);
    [mats,erra] = rla.get_fw_mats(kh,src_info,bctrans,sensor_info,opts);
    
    h_bd = src_info.h;
    src = zeros(4,n);
    src(1,:) = src_info.xs;
    src(2,:) = src_info.ys;
    src(3,:) = src_info.dxs;
    src(4,:) = src_info.dys;    

    S_tgt_out = slmat_out(zks(2,ik),h_bd,src,tgt(:,~in));
    D_tgt_out = dlmat_out(zks(2,ik),h_bd,src,tgt(:,~in)); 
    [m1,n1] = size(D_tgt_out);
    SD_tgt_out = zeros(m1,2*n1,'like',1.0+1i);
    SD_tgt_out(:,1:2:end) = -S_tgt_out/bs(2,ik);
    SD_tgt_out(:,2:2:end) = D_tgt_out/bs(2,ik);
    data_to_out_trans = SD_tgt_out*mats.inv_Fw_mat;    
    
    S_tgt_in = slmat_out(zks(1,ik),h_bd,src,tgt(:,in));
    D_tgt_in = dlmat_out(zks(1,ik),h_bd,src,tgt(:,in)); 
    [m1,n1] = size(D_tgt_in);
    SD_tgt_in = zeros(m1,2*n1,'like',1.0+1i);
    SD_tgt_in(:,1:2:end) = -S_tgt_in/bs(1,ik);
    SD_tgt_in(:,2:2:end) = D_tgt_in/bs(1,ik);
    data_to_in_trans = SD_tgt_in*mats.inv_Fw_mat;    


    % get bdry data and solutions
    
    t_dir = -pi/2;
    x_dir = cos(t_dir)';
    y_dir = sin(t_dir)';
    xs = src_info.xs(:)';
    ys = src_info.ys(:)';
    ds = src_info.ds(:)';
    dxs = src_info.dxs(:)';
    dys = src_info.dys(:)';
   
    fields = [];
    fields.uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    fields.dudninc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,1) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));

    q = 0.5*(as(1,ik)/bs(1,ik) + as(2,ik)/bs(2,ik));


    bd_data_trans = zeros(2*n,1);
    bd_data_trans(1:2:end,:) = -as(2,ik)*fields.uinc(:)/q;
    bd_data_trans(2:2:end,:) = -bs(2,ik)*fields.dudninc(:);          

    uscat_out_trans = data_to_out_trans*bd_data_trans;
    uscat_in_trans = data_to_in_trans*bd_data_trans;

    uinc_tgt = exp(1i *kh * (bsxfun(@times,tgt(1,:),x_dir) + ...
        bsxfun(@times,tgt(2,:),y_dir)));
    uinc_tgt = uinc_tgt(:);

    uplot_trans = zeros(size(xx)); 

    uplot_trans(~in) = uscat_out_trans + uinc_tgt(~in);
    uplot_trans(in) = uscat_in_trans;

    uplot_trans = reshape(uplot_trans,size(xx));

    fig = figure(1);
    fname = [fnamebase, sprintf('%d.png',kh)];
    clf
    h = pcolor(xx,yy,imag(uplot_trans)); set(h,'EdgeColor','none');
    clim([-2,2])
    colormap(brewermap([],'RdBu'));
    hold on
    plot(src_info.xs,src_info.ys,'k-','LineWidth',4);
    axis equal; axis tight;
    set(gca,"Visible","off")

    saveas(fig,fname);
    
end
