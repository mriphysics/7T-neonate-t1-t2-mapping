%%% 2024-10-01: Example T1 map script, working on phantom data
%%% Note that this data was acquired using a purposely lowered voltage for
%%% the inversion pulse, to create a situation with incomplete inversion
%%% Author: Shaihan Malik, King's College London

%% T1 map reconstruction
load data/t1data.mat
[nx ny ni] = size(imt1);

saveresults = 0;  % if 1 this runs optimization, if 0 it loads stored results

%% display the images 


figfp(2)
for ii=1:8
    subplot(2,4,ii)
    imsjm(imt1(:,:,ii))
    title(sprintf('Ti=%dms',ti(ii)),'fontsize',16)
end
colormap gray
spsjm('tgap',0.05,'gap',[0.01 0.05])
setpospap([50         332        1009         438])


%% Fitting of whole image  -- here estimate T1 with and without IE correction

%%% skip some echoes
echo_filter = [1 1 1 1 1 1 1 1]; 

mask = imt1(:,:,1)>300;

if saveresults
%%% Store all the results in a cell array. Array element 1 = with IE and
%%% cell array element 2 = without IE term

t1map={};
m0map={};
iemap={};
res1map={};

for IX = 1:2 % 1=with IE and 2=without IE

    %%% we set the IE term by setting the upper bound of the IE term in the
    %%% fit to either 0.5 (which means include it in the fit) or 0.0 (which
    %%% means do not include in the fit)

    switch IX
        case {1}
            eps_upper_bound = 0.5;
        case {2}
            eps_upper_bound = 0.0;
    end

    %%% perform fit
    t1map{IX} = zeros([nx ny]);
    m0map{IX} = zeros([nx ny]);
    iemap{IX} = zeros([nx ny]);
    res1map{IX}=zeros([nx ny]);

    opt = optimoptions('fmincon');
    opt.Display = "none";
    for ii=1:nx
        for jj=1:ny
            if mask(ii,jj)


                data = squeeze(imt1(ii,jj,:));


                sig = @(x)(x(1)*abs(1-2*(1-x(3))*exp(-ti(:)/x(2))));
                cf = @(x)(norm(echo_filter(:).*(data-sig(x)))^2);

                % guess at T1
                [~,mti]=min(data);
                t1guess = ti(mti)/log(2);
                xmin = fmincon(cf,[data(end) t1guess 0.1],[],[],[],[],[0 0 0],[data(end)*1.3 8000 eps_upper_bound],[],opt);

                m0map{IX}(ii,jj)=xmin(1);
                t1map{IX}(ii,jj)=xmin(2);
                iemap{IX}(ii,jj)=xmin(3);

                %%% residuals
                res1map{IX}(ii,jj) = cf(xmin);
            end
        end
        disp([ii IX])
    end

end

save outputs/t1map t1map m0map res1map iemap

else
    load outputs/t1map.mat
end


%% Now make image summarizing this data

figfp(10)

%%% get reverse gray colormap
cm = colormap('gray');
cm = flip(cm,1);

nr=2;nc=4;

fs=12;
fs1=12;
for IX=1:2

    subplot(nr,nc,1+(IX-1)*nc)
    imsjm(t1map{IX},[2000 3200],'rot',0)

    h2=imsjm(t1map{IX},[2000 3200]);
    h2.AlphaData=mask;

    colormap(gca,'inferno')
    axis off
    if IX==1
        title('T_1 (ms)','FontSize',fs1)
    else
        c1=colorbar;
        c1.Location="southoutside";
        c1.FontSize=fs;
        ylabel(c1,'ms')
    end

    subplot(nr,nc,2+(IX-1)*nc)
    h2=imsjm(abs(iemap{IX}),[0 0.15],'rot',0);
    h2.AlphaData=mask;
    colormap(gca,cm)
    axis off
    if IX==1
        title('Inversion inefficiency (\epsilon)','FontSize',fs1)
    else
        c2=colorbar;
        c2.Location="southoutside";
        c2.FontSize=fs;
        ylabel(c2,'au')
    end

    subplot(nr,nc,3+(IX-1)*nc)
    h2=imsjm(100*sqrt(res1map{IX}./sum(imt1.^2,3)),[0 10]);
    h2.AlphaData=mask;
    colormap(gca,'parula')
    axis off
    if IX==1
        title('residuals (%)','FontSize',fs1)
    else
        c3=colorbar;
        c3.Location="southoutside";
        c3.FontSize=fs;
        ylabel(c3,'%')
    end
end


subplot(nr,nc,[nc 2*nc])
hold on 
for ii=1:2
    histogram(t1map{ii}(mask))
end
grid on
xlabel('T_1 (ms)')
title('T_1 histograms')
legend('With \epsilon','Without \epsilon','location','northwest','fontsize',13)

gg = get(gcf,'children')
%%% 

% add text
axes(gg(11))
tt = text(-15,130,'With \epsilon','fontsize',fs1+3,'fontweight','bold','Rotation',90)
tt = text(-15,330,'Without \epsilon','fontsize',fs1+3,'fontweight','bold','Rotation',90)

gg(11).Position = [0.07 0.58 0.16 0.34];
gg(8).Position = [0.07 0.2 0.16 0.34];
gg(7).Position = [0.11 0.11 0.1 0.03];

xs=0.17;
gg(10).Position = [0.07 0.58 0.16 0.34]+[xs 0 0 0];
gg(6).Position = [0.07 0.2 0.16 0.34]+[xs 0 0 0];
gg(5).Position = [0.11 0.11 0.1 0.03]+[xs 0 0 0];

xs=0.34;
gg(9).Position = [0.07 0.58 0.16 0.34]+[xs 0 0 0];
gg(4).Position = [0.07 0.2 0.16 0.34]+[xs 0 0 0];
gg(3).Position = [0.11 0.11 0.1 0.03]+[xs 0 0 0];

gg(2).Position = [0.62 0.25 0.35 0.5];

setpospap([ 100 100 912 360])

print -dpng -r300 outputs/SuppInfo_S3.png
