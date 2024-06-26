%RECONSTRUCTIONS_LAMABV_WITHCONSTINITIAL
%
% try the starfish again with the constant model result
% as initial guess. here you'll have to change the output
% file names as appropriate in the defs of A and B below
%

clearvars

A = load('../imp-data/data-out/test_016_tensdata_impck_io_antbar3_phaseoff_25-Jun-2024 23:20:56.mat');
B = load('../imp-data/data-out/test_016_tensdata_impck.mat');

%%

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

linestyles = {'b--','b:','b-.'};

fac = 1.1;

fig = figure(1);
clf;
set(gcf,'Position',[0,0,800,600])

tiledlayout(1,2,"TileSpacing","tight");

st = B.src_info;
sf0 = A.inv_data_all0{end}.src_info_all{end};
omega0 = A.inv_data_all0{end}.kh;
sf1 = A.inv_data_all{end}.src_info_all{end};
omega1 = A.inv_data_all{end}.kh;
xt = st.xs; yt = st.ys;

xmax = max(xt); xmin = min(xt); xc = (xmax+xmin)/2;
dx = (xmax-xmin)/2; 
x1 = xc-dx*fac; x2 = xc+dx*fac;
ymax = max(yt); ymin = min(yt); yc = (ymax+ymin)/2;
dy = (ymax-ymin)/2; 
y1 = yc-dy*fac; y2 = yc+dy*fac;

xw = 0.1;
yw = (y2-y1)*xw/(x2-x1);
xa = -0.7-xw; xb = -0.7+xw; ya = -yw; yb = yw;

for j = 1:2
    t = nexttile;

    plot(xt,yt,'k-','LineWidth',2)
    hold on
        
    xs0 = sf0.xs; ys0 = sf0.ys;
    xs1 = sf1.xs; ys1 = sf1.ys;

    if (j == 1)
    plot(xs0,ys0,'g-.','LineWidth',2)
    plot(xs1,ys1,'b-.','LineWidth',2)
    rectangle('Position',[xa,ya,(xb-xa),(yb-ya)])
    set(gca,'XTick',[],'YTick',[]);
    h = gca;
    axis equal tight
        xlim([x1,x2]); ylim([y1,y2]);
        legend('obstacle', "constant model $\omega=" + num2str(omega0) +"$",...
            "curvature model $\omega=" + num2str(omega1) +"$");
        title('using initial guess')
    else
    plot(xs0,ys0,'gx','LineWidth',2,'MarkerSize',10)
    plot(xs1,ys1,'bx','LineWidth',2,'MarkerSize',10)
    set(gca,'XTick',[],'YTick',[]);
    h = gca;
    axis equal tight
    xlim([xa,xb]); ylim([ya,yb]);
        title('detail')
    end


    fontsize(gca,scale=1.5);

end


saveas(fig,'starfish_lamabv_with_const_initial.epsc');
