%%%%%%%%%%
%%% Create montage images of all maps
%%% Sort the data by age
%%%
%%% Shaihan Malik & Aiman Mahmoud, King's College London, 2024

%% Load in data structure
addpath('lib');
load data/cohort_data.mat

%% T1 maps

sub_n = (1:length(pmas))';

% Find indices of empty cells in t1dataie
t1_empty = cellfun(@isempty, t1dataie);

% Filter out the empty cells
t1data_filtered = t1dataie(~t1_empty);
sub_n = sub_n(~t1_empty);

% Initialize pma_t1 with pmas values and filter out rows corresponding to empty t1dataie cells
pma_t1 = [pmas(~t1_empty), zeros(length(pmas(~t1_empty)), 1)];

% Sort pma_t1 based on the first column in ascending order
[~, ixs] = sort(pma_t1(:,1), 'ascend');
pma_t1 = pma_t1(ixs, :);
sub_n = sub_n(ixs,:);

% Use the same sorting indices to sort t1data_filtered
t1data_sorted = t1data_filtered(ixs);

% Convert weeks to weeks and additional days
for i = 1:size(pma_t1, 1)
    currentWeeks = pma_t1(i,1);
    fullWeeks = floor(currentWeeks);
    additionalDays = round((currentWeeks - fullWeeks) * 7); % Correct conversion
    pma_t1(i, 2) = additionalDays; % Update the second column with additionalDays
    pma_t1(i, 1) = fullWeeks; % Update the first column with fullWeeks, if necessary
end

%%% filter and sort the IE data as well
% Filter out the empty cells
iedata_filtered = ie_maps(~t1_empty);
iedata_sorted = iedata_filtered(ixs); % sort to increasing age



%% Make a figure to show all data

figure(1);clf
% Define the number of rows and columns
numRows = 8;
numCols = 5;


for ii = 1:length(t1data_sorted)
   
    % mask
    m = t1data_sorted{ii}>0;

    % Create subplot with manual position
    subplot(numRows,numCols,ii);
    h1 = imagesc(t1data_sorted{ii}, [1500 4500]); % Assuming imsjm is similar to imshow
    h1.AlphaData = m;

    title({sprintf('#%d, PMA %d^{+%d}', sub_n(ii), pma_t1(ii,1), pma_t1(ii,2))});
end

% Add a single colorbar to the figure
cb = colorbar;
colormap('inferno')
yl=ylabel(cb, 'T_1 (ms)','Rotation',0);
yl.Position = [2.482 1120 0];
%caxis([min(z) max(z)])
set(cb, 'Position', [0.88 0.4 0.015 0.2],'FontSize',13);

spsjm('auto','rgap',0.15,'tgap',0.02,'gap',[0.0 0.02])

% setpospap([20         334        1072         686])

setpospap([20         1100        540     750])

print -dpng -r300 outputs/t1_data_montage.png

%% Now same again for epsilon

figfp(2);

% Define the number of rows and columns
numRows = 8;
numCols = 5;


for ii = 1:length(t1data_sorted)
   
    % mask
    m = t1data_sorted{ii}>0;

    % Create subplot with manual position
    subplot(numRows,numCols,ii);
    h1 = imagesc(iedata_sorted{ii}, [0 0.4]); % Assuming imsjm is similar to imshow
    h1.AlphaData = m;

     title({sprintf('#%d, PMA %d^{+%d}', sub_n(ii), pma_t1(ii,1), pma_t1(ii,2))});

end

% Add a single colorbar to the figure
cb = colorbar;
cmap = colormap('gray');
colormap(flip(cmap,1));
yl=ylabel(cb, '\epsilon','Rotation',0);
yl.Position = [1.247 -0.0533 0];
%caxis([min(z) max(z)])
set(cb, 'Position', [0.88 0.4 0.015 0.2],'FontSize',13);
yl.FontSize = 20;
spsjm('auto','rgap',0.15,'tgap',0.02,'gap',[0.0 0.02])

setpospap([20         1100        540     750])
print -dpng -r300 outputs/ie_data_montage.png

%% T2

sub_n = (1:length(pmas))';

% Find indices of empty cells in t1dataie
t2_empty = cellfun(@isempty, t2corrected);

% Filter out the empty cells
t2data_filtered = t2corrected(~t2_empty);
sub_n = sub_n(~t2_empty);

% Initialize pma_t1 with pmas values and filter out rows corresponding to empty t1dataie cells
pma_t2 = [pmas(~t2_empty), zeros(length(pmas(~t2_empty)), 1)];

% Sort pma_t1 based on the first column in ascending order
[~, ixs] = sort(pma_t2(:,1), 'ascend');
pma_t2 = pma_t2(ixs, :);
sub_n = sub_n(ixs,:);

% Use the same sorting indices to sort t1data_filtered
t2data_sorted = t2data_filtered(ixs);

% Convert weeks to weeks and additional days
for i = 1:size(pma_t2, 1)
    currentWeeks = pma_t2(i,1);
    fullWeeks = floor(currentWeeks);
    additionalDays = round((currentWeeks - fullWeeks) * 7); % Correct conversion
    pma_t2(i, 2) = additionalDays; % Update the second column with additionalDays
    pma_t2(i, 1) = fullWeeks; % Update the first column with fullWeeks, if necessary
end


figfp(2);

% Define the number of rows and columns
numRows = 7;
numCols = 5;

for ii = 1:length(t2data_sorted)
  
    mm = t2data_sorted{ii}>0;

    % Create subplot with manual position
    subplot(numRows,numCols,ii);
    h2 = imsjm(t2data_sorted{ii}, [20 180]); % Assuming imsjm is similar to imshow
    h2.AlphaData = mm;

    title({sprintf('#%d, PMA %d^{+%d}', sub_n(ii), pma_t2(ii,1), pma_t2(ii,2))});

end

cmapt2 = colorcet('L9');
colormap(cmapt2)
% Add a single colorbar to the figure
cb = colorbar;
yl = ylabel(cb, 'T_2 (ms)','Rotation',0);
yl.Position = [1.8444 10.4001 0];
set(cb, 'Position', [0.90 0.4 0.015 0.2],'FontSize',13);

spsjm('auto','rgap',0.1,'tgap',0.02,'gap',[0.0 0.02])

setpospap([20         1100        590         750])

print -dpng -r300 outputs/t2_data_montage.png


%% ROIS

rois_display = uneroded_rois;

sub_n = (1:length(pmas))';

% Find indices of empty cells in t1dataie
rois_empty = cellfun(@isempty, rois_display);
figure
% Filter out the empty cells
rois_filtered = rois_display(~rois_empty);
sub_n = sub_n(~rois_empty);

% Initialize pma_t1 with pmas values and filter out rows corresponding to empty t1dataie cells
pma_rois = [pmas(~rois_empty), zeros(length(pmas(~rois_empty)), 1)];

% Sort pma_t1 based on the first column in ascending order
[~, ixs] = sort(pma_rois(:,1), 'ascend');
pma_rois = pma_rois(ixs, :);
sub_n = sub_n(ixs,:);

% Use the same sorting indices to sort t1data_filtered
rois_sorted = rois_filtered(ixs);

% Convert weeks to weeks and additional days
for i = 1:size(pma_rois, 1)
    currentWeeks = pma_rois(i,1);
    fullWeeks = floor(currentWeeks);
    additionalDays = round((currentWeeks - fullWeeks) * 7); % Correct conversion
    pma_rois(i, 2) = additionalDays; % Update the second column with additionalDays
    pma_rois(i, 1) = fullWeeks; % Update the first column with fullWeeks, if necessary
end


figfp(3);

% Define the number of rows and columns
numRows = 8;
numCols = 5;

%%% custom color map
cmap = linspecer(8,'sequential');
% change some colours
cmap(4,:)=[0.7 0.4 0.6];

for ii = 1:length(rois_sorted)
   
    %%% manually remap the labels to match the figure in the main paper:
    tmp = zeros(size(rois_sorted{ii}));

    %%% indices to map 
    ix = [1 2 5 6 7 8 9 13];
    %%% 
    for jj=1:8
        tmp(rois_sorted{ii}==ix(jj))=jj;
    end
    
    % make a mask
    mm = tmp>0;

    % Create subplot with manual position
    subplot(numRows,numCols,ii);
    h3=imsjm(tmp,[1 8]); % Assuming imsjm is similar to imshow
   
    h3.AlphaData = mm;
    title({sprintf('#%d, PMA %d^{+%d}', sub_n(ii), pma_t1(ii,1), pma_t1(ii,2))});

  end

colormap(cmap);

spsjm('auto','rgap',0.05,'tgap',0.02,'gap',[0.0 0.02])

% setpospap([20         334        1072         686])

setpospap([20         334        490     750])


print -dpng -r300 outputs/segmentation_montage.png

