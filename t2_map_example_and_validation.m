%%% 2024-10-01: Example T2 map script, working on phantom data
%%% Phantom data acquired using 3D single SE and 2D TSE, reconstructed
%%% using exponential fitting and dictionary
%%% Author: Shaihan Malik, King's College London

saveresults = 0;  % if 1 this runs optimization, if 0 it loads stored results
addpath('lib');
load data/t2data.mat



%% Either run fit or load results


if saveresults
    
%%% first do 3D spin echo data
[n1 n2 ne] = size(im3d); % Note this is just the central slice from a 3D slab
t2fit3d = zeros([n1 n2]);
s0map3d = zeros([n1 n2]);


t2fit3d = zeros([n1 n2]);
s0map3d = zeros([n1 n2]);
res2map3d=zeros([n1 n2]);
opt = optimoptions('fmincon');
opt.Display = "none";

ww = diag([1 1 1 1 1]);
for ii=1:n1
    for jj=1:n2
        if mask3d(ii,jj)==0
            continue
        end
        
        data = double(squeeze(im3d(ii,jj,:)));
        
        %% nonlinear least squares
        sig = @(x)(x(1)*exp(-te3d(:)/x(2)));
        cf = @(x)(norm(ww*(data-sig(x)))^2);
        xmin = fmincon(cf,[data(1) 10],[],[],[],[],[0 0],[3*data(1) 500],[],opt);
        
        s0map3d(ii,jj)=xmin(1);
        t2fit3d(ii,jj)=xmin(2);
        
        %% residuals
        res2map3d(ii,jj) = cf(xmin);
        
        disp([ii jj]);
    end
end


%%% Now exponential fit on 2D maps.
t2fit2d = zeros([nx ny]);
s0map2d = zeros([nx ny]);
res2map2d=zeros([nx ny]);
opt = optimoptions('fmincon');
opt.Display = "none";

for ii=1:nx
    for jj=1:ny
        if mask2d(ii,jj)==0
            continue
        end

        data = double(squeeze(im2d(ii,jj,:)));

        %%% nonlinear least squares
        sig = @(x)(x(1)*exp(-te2d(:)/x(2)));
        cf = @(x)(norm(data-sig(x))^2);
        xmin = fmincon(cf,[data(1) 10],[],[],[],[],[0 0],[3*data(1) 500],[],opt);

        s0map2d(ii,jj)=xmin(1);
        t2fit2d(ii,jj)=xmin(2);

        %%% residuals
        res2map2d(ii,jj) = cf(xmin);
        
        disp([ii jj]);
    end
end

%%% save the results
save('outputs/t2fitresults.mat','t2fit2d','t2fit3d','s0map2d','s0map3d','res2map2d','res2map3d');
else

    %%% load in fits
   load outputs/t2fitresults.mat
end


%% Now for dictionary lookup we also need B1 map

% Load in B1 map.... this AFI was resampled using MIRTK command below
% mirtk transform-image SPH_1_AFI.nii.gz output_slice_.nii.gz -interp NN -target SPH_1_2D.nii.gz
nii = niftiread('data/b1map_t2data.nii.gz');
b1map = double(nii).' / 600;

%%% overlay b1 map mask
figfp(1)
imsjm(b1map)
colormap inferno
colorbar
hold on
vv=visboundaries(mask2d,'Color', 'white', 'LineWidth', 2, 'EnhanceVisibility', false);




%% 2024-09-16: Repeat above but use dictionary matching method

[nx ny] = size(t2fit2d);
t2mapd = zeros([nx ny]);
s0mapd = zeros([nx ny]);

%%% load in dictionary
load data/dict_20230629.mat

for ii=1:nx
    for jj=1:ny
        if mask2d(ii,jj)
            data = double(squeeze(im2d(ii,jj,:)));
            
            [t2tmp,s0tmp] = dict_match(data,b1map(ii,jj),dict);

            s0mapd(ii,jj)=s0tmp;
            t2mapd(ii,jj)=t2tmp;
        end
    end

end

%% New figure with dictionary match result


%%% pad the 3D result with zeroes so it looks the same size as the 3D
%%% result
pp=6;
t2fit3dp = padarray(t2fit3d,[pp pp]);
mask3dp = padarray(mask3d,[pp pp]);

w1 = [40 200];
w2 = [40 200];


nr = 2;
nc = 3;
fs = 11;
fs2 = 11;

figfp(2)
subplot(nr,nc,1)
h1 = imsjm(t2fit3dp.*mask3dp,w1);
h1.AlphaData=mask3dp;

title('T_2 3D Spin Echo','fontsize',fs)
axis off

subplot(nr,nc,2)
h2 = imsjm(t2fit2d,w2);
h2.AlphaData=mask2d;

title('T_2 2D exponential fit','fontsize',fs)
axis off


subplot(nr,nc,3)
h3 = imsjm(t2mapd,w1);
h3.AlphaData=mask2d;

title('T_2 dictionary recon','fontsize',fs)
axis off

cc = colorbar;
cc.FontSize = fs;
cmap = colorcet('L9');
colormap(cmap)

%%% ROI analysis
subplot(nr,nc,4:6)
hold on
nbin=64;
w3 = [40 200];
histogram(t2fit2d(mask2d),nbin,'BinLimits',w3,'Normalization','probability');
histogram(t2mapd(mask2d),nbin,'BinLimits',w3,'Normalization','probability');
histogram(t2fit3d(mask3d),nbin,'BinLimits',w3,'Normalization','probability');
legend('2D exponential fit','2D dictionary recon','3D spin echo (reference method)')
grid on
xlabel('T_2 measurement')
ylabel('fraction of pixels')
set(gca,'FontSize',fs2)

gg = get(gcf,'Children');
ww=700;
hh=450;
setpospap([200 200 ww hh])

aspectratio = nx/ny;
hhh=0.42;

% hhh*hh/(www*ww) = aspectratio
www = (hhh*hh)/aspectratio/ww;
gg(6).Position = [0.01 0.5 www hhh];
gg(5).Position = [0.31 0.5 www hhh];
gg(4).Position = [0.61 0.5 www hhh];
gg(2).Position = gg(2).Position + [0 0.02 0 0];

% print -dpng -r300 phantom_validation_figure_dictionary_recon.png



%% print out stats:
fprintf(1,'3D SE: \t\tT2 = %2.1f\t±\t%2.1f ms\n',mean(t2fit3d(mask3d)),std(t2fit3d(mask3d)));
fprintf(1,'2D exp.: \tT2 = %2.1f\t±\t%2.1f ms\n',mean(t2fit2d(mask2d)),std(t2fit2d(mask2d)));
fprintf(1,'2D dict.: \tT2 = %2.1f\t±\t%2.1f ms\n',mean(t2mapd(mask2d)),std(t2mapd(mask2d)));

