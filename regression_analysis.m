%%% ROI analysis for infant data. Load in ROI measurements for each infant,
%%% over different tissue labels. Note that the data contains 61 infant
%%% datasets but not all contain T1 and T2 mapping (so are empty) and it
%%% also contains 13 tissue labels, though we only consider 8 of these.
%%% 2024-10-08: Shaihan Malik & Aiman Mahmoud, King's College London
%%% 2024-10-24: Update to use a mixed effects model

addpath('lib');

%%% Load in mat file containing the data
load data/cohort_data.mat

%%% Define labels for the ROIs we are interested in
all_labs = {'Cortical GM','White Matter','Lateral Ventricle',...
    'Cavum septum pellucidum','Brainstem','Cerebellum','Cerebellar Vermis','Basal Ganglia','Thalamus','Third Ventricle','Fourth Ventricle', 'CSF', 'PVFWM'};
%%% These ROIs are to be displayed from the 13 labels included in segmentation
roi_idx = [1,2,5,6,7,8,9,13];
newlabs = all_labs(roi_idx);

%%% Go through the T1 and T2 maps and save the average figures within each
%%% of the ROIs indexed above
N = length(t1dataie); % 61 total subjects

t1m_all = zeros([N 8]); % N subjects, 8 labels
t2m_all = zeros([N 8]); 
%%% fill the array with NaN - this is a label for missing data later
t1m_all(:) = NaN;t2m_all(:) = NaN;
for ii=1:N
    for jj=1:8
        ix = find(new_rois_t1{ii} == roi_idx(jj));
        tmp_t1 = t1dataie{ii}(ix);
        % remove NaNs before taking averages
        tmp_t1(isnan(tmp_t1)) = [];
        t1m_all(ii,jj) = median(tmp_t1);
        
        % now look at T2
        if ~isempty(t2corrected{ii})
            ix = find(new_rois_t2{ii} == roi_idx(jj));
            tmp_t2 = t2corrected{ii}(ix);
            tmp_t2(isnan(tmp_t2)) = [];
            t2m_all(ii,jj) = median(tmp_t2);
        end
    end
end


%%
figfp(1)

%%% Build samples for plotting. T1 first.
t1m = t1m_all;
% These scans are not included in regression
t1_rejection = [1 3 5 33];
t2_rejection = [1 3 5 33 49]; % additionally exclude 49 because of severe motion artefacts

%%% Remove the rejected subjects from variable used in linear regression
%%% (mark as NaN)
t1m(t1_rejection,:) = NaN;

%%% Create another variable for plotting of excluded scans
t1m_outliers = t1m_all(t1_rejection,:);
pmas_outliers_t1 = pmas(t1_rejection);
pnas_outliers_t1 = pmas(t1_rejection)-gabs(t1_rejection);


%%% Now the same thing but for T2
t2m = t2m_all;

%%% Remove the rejected ones from linear regression (set as Nan to flag)
t2m(t2_rejection,:) = NaN; 

%%% Create another variable for plotting of these
t2cm_outliers = t2m_all(t2_rejection,:);
pmas_outliers_t2 = pmas(t2_rejection);
pnas_outliers_t2 = pmas(t2_rejection)-gabs(t2_rejection);

%%% Indices of repeat scans
idx_rep = [[7 18];[13 14];[20 27];[21 28];[25 26];[50 51]];

%%% create unique patient ID
patid = 1:length(pmas);
patid(idx_rep(:,2)) = idx_rep(:,1);

%%% variable to save data 
plot_data = {};

%%% ===================================
%%% Perform regressions and create scatter plots all within this loop.

%%% do two regression analyses - single and multiple linear regression
lm_single = {};
lm_multiple = {}; %<-- for linear models including PNA

%%% params
mrkrsz = 30; % scatter plot marker size
%%% colormap for plotting - use COLORCET function    cmap = colorcet('L6');
cmap = colorcet('L8');
%%% stretch non-linear
x0 = linspace(0,1,256);
cmap = interp1(x0(:),cmap,x0(:).^0.7);

% range for PNA colours
pna_min = 0;
pna_max = 14;

% colour for outlier points
out_color = [0 0.7 0];

% Variable to save marginal R^2
Rsq_marg = {};

for jj=1:8 % loop over ROIs
     
    %%% perform linear model fit fot T1
    x = pmas-40; %%%<--- ref to 40wk (i.e. PMA-40 so intercept is value at 40wk)
    y = t1m(:,jj);
    z = pmas - gabs; %%% 3rd variable is PNA 
    id = patid(:);

    %%% for outliers plot
    k = t1m_outliers(:,jj);
    x2 = pmas_outliers_t1-40;
    z2 = pnas_outliers_t1;

    %find NaN values and remove them from x and y
    % Remove NaN values
    nanIndices={};
    nanIndices = isnan(y);
    x(nanIndices)=[];
    y(nanIndices)=[];
    z(nanIndices)=[];
    id(nanIndices) = [];


    %%% Put data into a table
    data_tab = table(x,y,z,id);
    data_tab.Properties.VariableNames = {'PMA','T1','PNA','sub_id'};

    %%% now use fitlme to fit mixed effects models. Do this either for
    %%% single or multiple regression
    lm_single{jj,1} = fitlme(data_tab, "T1 ~ PMA + (1|sub_id)", "FitMethod", "REML");

    % multiple regression also considers PNA
    lm_multiple{jj,1} = fitlme(data_tab, "T1 ~ PMA + PNA + (1|sub_id)", "FitMethod", "REML");

    % x range for confidence intervals
    xv = linspace(-7,13, 150);

    % Need to create a dummy table for the 'predict' function
    tblnew = table();
    tblnew.PMA = xv';
    tblnew.PNA = linspace(min(z),max(z), 150)';
    tblnew.sub_id = ones([150 1]);
    % now predict the conf intervals and plot value
    [ypred,yCI,DF] = predict(lm_single{jj,1},tblnew);

    %%% Colours for scatter plot
    cols = interp1(linspace(pna_min,pna_max,256), cmap, z, 'linear');
    cols2 = interp1(linspace(pna_min,pna_max,256), cmap, z2, 'linear');
    %%% values larger than max get assigned max color value
    iMax = find(z>pna_max);
    cols(iMax,:) = repmat(cmap(end,:),[length(iMax) 1]);
    iMax = find(z2>pna_max);
    cols2(iMax,:) = repmat(cmap(end,:),[length(iMax) 1]);
    

    subplot(4,4,jj)
    hold on
    patch([xv fliplr(xv)], [yCI(:,1);flip(yCI(:,2),1)], zeros(300,1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5, 'EdgeColor','none')
    p2{jj} = plot(xv,ypred);
    ss = scatter(x, y,mrkrsz,cols,'filled');
    set(ss,'MarkerEdgeColor',[0 0 0],'LineWidth',0.5);

    ss2 = scatter(x2,k,mrkrsz,cols2,'filled','^');
    set(ss2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',out_color,'LineWidth',0.5);

    hold off
    grid on
    xlim([34 54]-40)
    ylim([1600 3800])
    if jj>4
        xlabel('PMA (weeks)');% only x-axis for lower plots
    else
        xlabel('');
    end
    ylabel('T_1 (ms)')
    legend('off')
    title(sprintf('%s T_1',newlabs{jj}),'fontsize',11);%
    set(gca,'XTick',-5:5:10,'XTickLabel',{'35','40','45','50'}); % label as true PMA (not minus 40)

    p2{jj}.LineWidth = 1.;
    p2{jj}.Color = [0 0. 0.];
    ss.MarkerEdgeColor=[0 0 0];

    %%% Compute the marginal Rsquared
    Rsq_marg{jj,1} = 1-var(residuals(lm_single{jj,1},'Condition',0))/var(data_tab.T1);

    % add text
    if jj<8
        %tt = text(40-40,3500,sprintf('N=%d;\tR^2=%1.2f',lm_single{jj,1}.NumObservations,Rsq_marg{jj,1}),'color',[0 0 1]);
        tt = text(45-40,3500,sprintf('N=%d',lm_single{jj,1}.NumObservations),'color',[0 0 1]);
    else
        %tt = text(40-40,2000,sprintf('N=%d;\tR^2=%1.2f',lm_single{jj,1}.NumObservations,Rsq_marg{jj,1}),'color',[0 0 1]);
        tt = text(45-40,2000,sprintf('N=%d',lm_single{jj,1}.NumObservations),'color',[0 0 1]);
    end
  
    %%% NOW THE SAME FOR T2

    %%% perform linear model fit
    x = pmas-40; %%% <--- reference to 40 wks
    y = t2m(:,jj);
    z = pmas-gabs;
    id = patid(:);

    % Remove NaN values
    nanIndices={};
    nanIndices = isnan(y);
    x(nanIndices)=[];
    y(nanIndices)=[];
    z(nanIndices)=[];
    id(nanIndices)=[];

    %%% for outliers plot
    k = t2cm_outliers(:,jj);
    x2 = pmas_outliers_t2-40;
    z2 = pnas_outliers_t2;

    %%% Put data into a table
    data_tab = table(x,y,z,id);
    data_tab.Properties.VariableNames = {'PMA','T2','PNA','sub_id'};

    %%% now use fitlme to fit mixed effects models. Do this either for
    %%% single or multiple regression
    lm_single{jj,2} = fitlme(data_tab, "T2 ~ PMA + (1|sub_id)", "FitMethod", "REML");

    % multiple regression also considers PNA
    lm_multiple{jj,2} = fitlme(data_tab, "T2 ~ PMA + PNA + (1|sub_id)", "FitMethod", "REML");

    % x range for confidence intervals
    xv = linspace(-8,5, 150);

    % Need to create a dummy table for the 'predict' function
    tblnew = table();
    tblnew.PMA = xv';
    tblnew.PNA = linspace(min(z),max(z), 150)';
    tblnew.sub_id = ones([150 1]);
    % now predict the conf intervals and plot value
    [ypred,yCI,DF] = predict(lm_single{jj,2},tblnew);

    %%% interpolate colors
    cols = interp1(linspace(pna_min, pna_max,256), cmap, z, 'linear');
    cols2 = interp1(linspace(pna_min,pna_max,256), cmap, z2, 'linear');
    %%% values larger than max get assigned max color value
    iMax = find(z>pna_max);
    cols(iMax,:) = repmat(cmap(end,:),[length(iMax) 1]);
    iMax = find(z2>pna_max);
    cols2(iMax,:) = repmat(cmap(end,:),[length(iMax) 1]);

    subplot(4,4,jj+8)
    hold on
    patch([xv fliplr(xv)], [yCI(:,1);flip(yCI(:,2),1)], zeros(300,1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5, 'EdgeColor','none')
    p2{jj} = plot(xv,ypred);
    ss = scatter(x, y,mrkrsz,cols,'filled');
    set(ss,'MarkerEdgeColor',[0 0 0],'LineWidth',0.5);

    ss2 = scatter(x2,k,mrkrsz,cols2,'filled','^');
    set(ss2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',out_color,'LineWidth',0.5);

    hold off
    grid on
    xlim([32 45]-40)

    %%% Compute the marginal Rsquared
    Rsq_marg{jj,2} = 1-var(residuals(lm_single{jj,2},'Condition',0))/var(data_tab.T2);


    %%% use two different y-limits
    switch jj
        case {1,2,4,8}
            ylim([70 200]);
            %tt = text(35-40,190,sprintf('N=%d;\tR^2=%1.2f',lm_single{jj,2}.NumObservations,Rsq_marg{jj,2}),'color',[0 0 1]);
            tt = text(42-40,190,sprintf('N=%d',lm_single{jj,2}.NumObservations),'color',[0 0 1]);
        case {3,5,6,7}
            ylim([70 120]);
            %tt = text(35-40,115,sprintf('N=%d;\tR^2=%1.2f',lm_single{jj,2}.NumObservations,Rsq_marg{jj,2}),'color',[0 0 1]);
            tt = text(42-40,115,sprintf('N=%d',lm_single{jj,2}.NumObservations),'color',[0 0 1]);
    end
    if jj>4
        xlabel('PMA (weeks)');% only x-axis for lower plots
    else
        xlabel('');
    end
    ylabel('T_2 (ms)')
    title(sprintf('%s T_2',newlabs{jj}),'fontsize',11)

    set(gca,'XTick',-5:5:10,'XTickLabel',{'35','40','45','50'});
    %%% remove legend
    legend('off')
    p2{jj}.LineWidth = 1.;
    p2{jj}.Color = [0 0. 0.];
    ss.MarkerEdgeColor=[0 0 0];

end

setpospap([163 1050 800 750])
%%%
gg = getAxesChildren(gcf);

%%% move the T2 plots down a bit and left
for ii=1:2:16
    pos = gg(ii).Position;
    pos(2) = pos(2) - 0.05;
    pos(1) = pos(1) - 0.04;
    gg(ii).Position = pos;
end

%%% move the T1 plots up a bit and left
for ii=2:2:16
    pos = gg(ii).Position;
    pos(2) = pos(2) + 0.03;
    pos(1) = pos(1) - 0.04;
    gg(ii).Position = pos;
end

% Add a single colorbar to the figure
cb = colorbar;
colormap(cmap);
ylabel(cb, 'Postnatal Age (PNA, weeks)');
caxis([pna_min pna_max])
set(cb, 'Position', [0.91 0.3 0.015 0.4], 'Units', 'normalized');
cb.FontSize = 12;


%%% Add some annotations for parts A and B
axes(gg(8))
text(-15,1650,'(a)','FontSize',24,'FontWeight','bold')
axes(gg(7))
text(-14,68,'(b)','FontSize',24,'FontWeight','bold')


print -dpng -r300 outputs/t1_t2_combined.png

%% compile the result of single linear regression


%%% put into a cell array
t1t2_regression_single = {};

for ii=1:8 % ROIs
    for jj=1:2 % T1 or T2
        
        %%% T1 values
        cci = lm_single{ii,jj}.coefCI;
        cci = round(cci);
        t1t2_regression_single{ii,1,jj} = sprintf('%d (%d,%d)',round(lm_single{ii,jj}.Coefficients.Estimate(1)),cci(1,1),cci(1,2));
            
        % coeff for PMA
        switch jj
            case 1 % T1 - round to integer
                t1t2_regression_single{ii,2,jj} = sprintf('%d (%d,%d)',round(lm_single{ii,jj}.Coefficients.Estimate(2)),cci(2,1),cci(2,2));
            case 2 % T2 - round to 1dp
                t1t2_regression_single{ii,2,jj} = sprintf('%1.1f (%1.1f,%1.1f)',(lm_single{ii,jj}.Coefficients.Estimate(2)),cci(2,1),cci(2,2));
        end
        
        
        % p-val
        pval = lm_single{ii,jj}.Coefficients.pValue(2);
        if pval < 0.001
            t1t2_regression_single{ii,3,jj} = '< 0.001';
        else
            t1t2_regression_single{ii,3,jj} = sprintf('%1.3f',pval);
        end

        %  Rsquared
        t1t2_regression_single{ii,4,jj} = sprintf('%1.2f',Rsq_marg{ii,jj});

    end
end

writecell(t1t2_regression_single(:,:,1),'outputs/regression_single.xlsx','Sheet','T1')
writecell(t1t2_regression_single(:,:,2),'outputs/regression_single.xlsx','Sheet','T2')


%% Put the result of multiple linear regression into an excel spreadsheet

%%% put into a cell array
t1t2_regression = {};

for ii=1:8 % ROIs
    for jj=1:2 % T1 or T2
        
        %%% intercept value
        cci = lm_multiple{ii,jj}.coefCI;
        cci = round(cci);
        t1t2_regression{ii,1,jj} = sprintf('%d (%d,%d)',round(lm_multiple{ii,jj}.Coefficients.Estimate(1)),cci(1,1),cci(1,2));
   
        % coeff for PMA
        switch jj
            case 1 % T1 - round to integer
                t1t2_regression{ii,2,jj} = sprintf('%d (%d,%d)',round(lm_multiple{ii,jj}.Coefficients.Estimate(2)),cci(2,1),cci(2,2));
            case 2 % T2 - round to 1dp
                t1t2_regression{ii,2,jj} = sprintf('%1.1f (%1.1f,%1.1f)',(lm_multiple{ii,jj}.Coefficients.Estimate(2)),cci(2,1),cci(2,2));
        end

        % p-val
        pval = lm_multiple{ii,jj}.Coefficients.pValue(2);
        if pval < 0.001
            t1t2_regression{ii,3,jj} = '< 0.001';
        else
            t1t2_regression{ii,3,jj} = sprintf('%1.3f',pval);
        end

        % coeff for PMA
        switch jj
            case 1 % T1 - round to integer
                t1t2_regression{ii,4,jj} = sprintf('%d (%d,%d)',round(lm_multiple{ii,jj}.Coefficients.Estimate(3)),cci(3,1),cci(3,2));
            case 2 % T2 - round to 1dp
                t1t2_regression{ii,4,jj} = sprintf('%1.1f (%1.1f,%1.1f)',(lm_multiple{ii,jj}.Coefficients.Estimate(3)),cci(3,1),cci(3,2));
        end
        
        
        
        % p-val
        pval = lm_multiple{ii,jj}.Coefficients.pValue(3);
        if pval < 0.001
            t1t2_regression{ii,5,jj} = '< 0.001';
        else
            t1t2_regression{ii,5,jj} = sprintf('%1.3f',pval);
        end

        % Compute the marginal Rsquared for this
        t1t2_regression_single{ii,4,jj} = 1-var(residuals(lm_multiple{ii,jj},'Condition',0))/var(lm_multiple{ii,jj}.response);
    end
end

writecell(t1t2_regression(:,:,1),'outputs/regression_multiple.xlsx','Sheet','T1')
writecell(t1t2_regression(:,:,2),'outputs/regression_multiple.xlsx','Sheet','T2')


