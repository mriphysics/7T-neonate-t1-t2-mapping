%%% Example for T1 and T2 estimation plus display of maps, for one example
%%% neonatal dataset

addpath('lib');

%%% data has been pre-saved into .mat to ensure anonymized
load data/neonate_data_example.mat

%%% Flag to either compute the maps, or load from data
saveresults = 1; % 0 - don't calculate, 1 - calculate and save

%%% mask for fitting
mask = imt2(:,:,1)>100;

%% display the images 

figfp(1)
for ii=1:3
    subplot(1,3,ii)
    imsjm(imt2(:,:,ii),[0 1500])
    title(sprintf('TE=%dms',te(ii)),'fontsize',20)
end
colormap gray
spsjm('tgap',0.1)
setpospap([50 100 866 340])

figfp(2)
for ii=1:8
    subplot(2,4,ii)
    %imsjm(imt1(:,:,ii),[0 1500])
    imsjm(imt1(:,:,ii))
    title(sprintf('Ti=%dms',ti(ii)),'fontsize',16)
end
colormap gray
spsjm('tgap',0.05,'gap',[0.01 0.05])
setpospap([50 100 900 600])


%% T1 estimation

if saveresults

%%% perform fit
t1map = zeros([nx ny]);
m0map = zeros([nx ny]);
iemap = zeros([nx ny]);
res1map=zeros([nx ny]);

opt = optimoptions('fmincon');
opt.Display = "none";
for ii=1:nx
    for jj=1:ny
        if mask(ii,jj)

            data = squeeze(imt1(ii,jj,:));
            sig = @(x)(x(1)*abs(1-2*(1-x(3))*exp(-ti(:)/x(2))));
            cf = @(x)(norm((data-sig(x)))^2);

            % guess at T1
            [~,mti]=min(data);
            t1guess = ti(mti)/log(2);
            xmin = fmincon(cf,[data(end) t1guess 0.],[],[],[],[],[0 0 0],[data(end)*1.3 6000 0.2],[],opt);

            m0map(ii,jj)=xmin(1);
            t1map(ii,jj)=xmin(2);
            iemap(ii,jj)=xmin(3);

            %%% residuals
            res1map(ii,jj) = cf(xmin);
        end
    end
    disp([ii])
end

save data/neonate_t1map.mat t1map m0map res1map iemap

else
    load data/neonate_t1map
end


%% T2 dictionary recon

if saveresults

    %%% Simple dictionary lookup
    load data/dict_20230629.mat
    %%% perform fit

    t2mapd = zeros([nx ny]);
    s0mapd = zeros([nx ny]);

    for ii=1:nx
        for jj=1:ny
            if mask(ii,jj)

                data = squeeze(imt2(ii,jj,:));

                [t2tmp,s0tmp] = dict_match(data,b1map(ii,jj),dict);

                s0mapd(ii,jj)=s0tmp;
                t2mapd(ii,jj)=t2tmp;
            end
        end

    end

    %%% For this one only, also do dictionary lookup including T1
    load data/dict_incT1_20230629.mat
    %%% perform fit

    t2mapdt1 = zeros([nx ny]);
    s0mapdt1 = zeros([nx ny]);

    for ii=1:nx
        for jj=1:ny
            if mask(ii,jj)

                data = squeeze(imt2(ii,jj,:));

                [t2tmp,s0tmp] = dict_match_inct1(data,b1map(ii,jj),t1map(ii,jj),dict);

                s0mapdt1(ii,jj)=s0tmp;
                t2mapdt1(ii,jj)=t2tmp;
            end
        end

    end

    %%% save all data
    save data/neonate_t2map.mat t2mapd s0mapd t2mapdt1 s0mapdt1

else
    load data/neonate_t2map.mat
end



%% Voxel-wise fitting

%%% skip some echoes
echo_filter = [1 1 1 1 1 1 1 1];

% Display one of the images, and then allow the user to click on one pixel
figure(10);clf
set(gcf,'Visible','on')
subplot(2,1,1)
imagesc(sum(imt1,3));colormap gray
title('Click to select a voxel')

%% Now make some images

figfp(10)

nr=2;nc=2;

subplot(nr,nc,1)
imsjm(m0map,'rot',0)
colormap(gca,'gray')
title('S_0 image','FontSize',20)

subplot(nr,nc,2)
imsjm(t1map,[1500 4000],'rot',0)
colorbar
colormap(gca,'inferno')
c1=colorbar;
c1.Location="southoutside";
title('T_1 (ms)','FontSize',20)
c1.Position = [0.63 0.52 0.2500 0.0200];
c1.FontSize=13;

subplot(nr,nc,3)
imsjm(abs(iemap),[0 0.3],'rot',0)
colormap(gca,'parula')
c2=colorbar;
c2.Location="southoutside";
title('\epsilon','FontSize',30)
c2.Position = [0.12 0.05 0.25 0.02];
c2.FontSize=13;

subplot(nr,nc,4)
imsjm(100*sqrt(res1map./sum(imt1.^2,3)),[0 15])
colormap(gca,'gray')
title('residuals','FontSize',20)
c3=colorbar;
c3.Location="southoutside";
title('residuals (%)','FontSize',20)
c3.Position = [0.63 0.05 0.25 0.02];
c3.FontSize=13;


spsjm('bgap',0.1,'tgap',0.1,'gap',[0.01 0.1])

setpospap([100 1 500 796])

% print -dpng -r300 T1map.png



%% T1 and T2 maps combined

figfp(1)

nr=1;nc=2;

subplot(nr,nc,1)
imsjm(t1map,[1600 4000],'rot',0)
colorbar
colormap(gca,'inferno')
c1=colorbar;
c1.Location="southoutside";
title('T_1 (ms)','FontSize',20)
c1.Position = [0.12 0.05 0.25 0.03];
c1.FontSize=13;

subplot(nr,nc,2)
imsjm(t2map,[50 270],'rot',0)
colorbar
colormap(gca,'viridis')
c1=colorbar;
c1.Location="southoutside";
title('T_2 (ms)','FontSize',20)
c1.Position = [0.62 0.05 0.25 0.03];
c1.FontSize=13;



spsjm('bgap',0.1,'tgap',0.1)

setpospap([100 100 650 400])

print -dpng -r300 t1_t2_neo.png


%% Include a T1 corrected dictionary lookup
%%% Get B1
nii = niftiread('afi_resampled.nii.gz');
b1map = double(nii) / 600;
b1map = imrotate(b1map,90); 

%%% get dict
load dict_incT1_20230629.mat
load t1map_sept24.mat
%%% perform fit

t2mapdt1 = zeros([nx ny]);
s0mapdt1 = zeros([nx ny]);

for ii=1:nx
    for jj=1:ny
        if mask(ii,jj)

            data = squeeze(imt2(ii,jj,:));
            
            [t2tmp,s0tmp] = dict_match_inct1(data,b1map(ii,jj),t1map(ii,jj),dict);

            s0mapdt1(ii,jj)=s0tmp;
            t2mapdt1(ii,jj)=t2tmp;
        end
    end

end

%% T1 map figure

%%% zoom in a bit
xr = 25:165;
yr = 20:155;

figure(2);clf;
nr=1;nc=3;
fs=16;

subplot(nr,nc,1)
h2=imsjm(t1map(xr,yr),[1500 4500]);
h2.AlphaData=mask(xr,yr);
colorbar
% cm = colorcet('L8');
colormap(gca,'inferno')
c1=colorbar;
c1.Location="eastoutside";
title('T_1','FontSize',fs)
c1.Position = [0.275 0.3 0.02 0.44];
c1.FontSize=13;
ylabel(c1,'ms','Position',[0.5 1300],'Rotation',0,'FontSize',16,'FontWeight','bold')

subplot(nr,nc,2)
h3=imsjm((iemap(xr,yr)),[0 0.4]);
h3.AlphaData=mask(xr,yr);
cm = colorcet('L1');
cm = flip(cm,1);
colormap(gca,cm);
c2=colorbar;
c2.Location="eastoutside";
title('\epsilon','FontSize',30)
c2.Position = [0.60 0.3 0.02 0.44];
c2.FontSize=13;

subplot(nr,nc,3)
h6=imsjm(b1map(xr,yr),[0.4 1.2]);
h6.AlphaData=mask(xr,yr);
colormap(gca,parula)
cc6=colorbar;
hold on
title('B_1^{rel}','fontsize',fs,'fontangle','italic')

cc6.Position = [0.925 0.3 0.02 0.44];
cc6.FontSize=13;

setpospap([100 100 935 300]);

spsjm('auto','tgap',0.1,'gap',[0.05 0.1],'rgap',0.08)

aa = text(-310,130,'(a)','fontsize',26,'FontWeight','bold');
bb = text(-150,130,'(b)','fontsize',26,'FontWeight','bold');
cc = text(10,130,'(c)','fontsize',26,'FontWeight','bold');

print -dpng -r300 outputs/T1map_neonate_figure.png


%% T2 map figure comparing with and without T1

win = [20 180];
win2 = [-3 3];
figure(2)
clf

%%% zoom in a bit
xr = 25:165;
yr = 20:155;

%%% colormaps
cmapd1 = colorcet('D1');
cmapt2 = colorcet('L9');%colormap(viridis);%
cmapratio = colorcet('L4');
fs = 15;


subplot(131)
h2=imsjm(t2mapd(xr,yr),win);
h2.AlphaData=mask(xr,yr)&~isnan(t2mapd(xr,yr));
title('T_2 (fixed T_1)','fontsize',fs)
cc1=colorbar;
colormap(gca,cmapt2)

subplot(132)
h3=imsjm(t2mapdt1(xr,yr),win);
h3.AlphaData=mask(xr,yr)&~isnan(t2mapdt1(xr,yr));
title('T_2 (incl. T_1 map)','fontsize',fs)
cc2=colorbar;
colormap(gca,cmapt2)


subplot(133)
h5=imsjm(t2mapd(xr,yr)-t2mapdt1(xr,yr),win2);
h5.AlphaData=mask(xr,yr)&~isnan(t2mapdt1(xr,yr));
colormap(gca,cmapd1);
cc3=colorbar;
title('\DeltaT_2','fontsize',fs)

setpospap([100 100 960 310]);

spsjm('auto','tgap',0.1,'gap',[0.05 0.1],'rgap',0.08)

%
cc1.Position = [0.27 0.3 0.015 0.4];
cc1.FontSize = 12;
ylabel(cc1,'ms','Position',[0.5 15],'Rotation',0,'FontSize',14,'FontWeight','bold')

cc2.Position = [0.6 0.3 0.015 0.4];
cc2.FontSize = 12;
ylabel(cc2,'ms','Position',[0.5 15],'Rotation',0,'FontSize',14,'FontWeight','bold')

cc3.Position = [0.93 0.3 0.015 0.4];
cc3.FontSize = 12;
ylabel(cc3,'ms','Position',[0.5 -3],'Rotation',0,'FontSize',14,'FontWeight','bold')



%
aa = text(-310,125,'(a)','fontsize',26,'FontWeight','bold');
ab = text(-150,125,'(b)','fontsize',26,'FontWeight','bold');
ac = text(10,125,'(c)','fontsize',26,'FontWeight','bold');


print -dpng -r300 outputs/T2map_neonate_figure.png


