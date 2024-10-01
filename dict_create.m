%%% 2024-09-16: construct a dictionary for reconstruction of the T2 data.
%%% This is based on a Bloch simulation performed separately (code not
%%% included for that; details in manuscript).
%%%
%%% This script creates dictionaries and evaluates the effect of including
%%% and excluding T1 in the lookup.
%%% Shaihan Malik, King's College London, 2024


%%% Load in raw simulated data
% S is a 4D matrix with dimensions nb1 x nt2 x nt1 x nechoes
load data/S_isochromat_20230629.mat

% Now, extract only the echo times measured 
teidx = [5 13 24]; %<--- these are echoes to use
ne = 3;

% take absolute signal at each echo time
Sdict = abs(S(:,:,:,teidx));

%% First interpolate the dictionary to T1=2.6s because we will fix there

Sdict_fixT1 = zeros([nb1 nt2 ne]);

for ii=1:ne
    Sdict_fixT1(:,:,ii) = interp3(T2,b1fact,T1, ...
        Sdict(:,:,:,ii),T2,b1fact,2600);
end

[t2mesh,b1mesh] = meshgrid(T2,b1fact);

%% Up-sample this dictionary to make more values using linear interpolation

%%% upsample to 1ms and 0.01 B1rel accuracy
T2hi = min(T2(:)):1:max(T2(:));
b1hi = min(b1fact(:)):0.01:max(b1fact(:));
nt2hi = length(T2hi);
nb1hi = length(b1hi);
[t2mesh_hi,b1mesh_hi] = meshgrid(T2hi,b1hi);

%%% Upsample dictionary
Sdict_hi = zeros([nb1hi nt2hi ne]);
for ii=1:ne
    Sdict_hi(:,:,ii) = interp2(t2mesh,b1mesh,Sdict_fixT1(:,:,ii),t2mesh_hi,b1mesh_hi);
end

%%% normalise the dictionary through the echoes
Snorm = vecnorm(Sdict_hi,2,3); % 2-norm in 3rd direction
Sdict_hi_norm = Sdict_hi ./ repmat(Snorm,[1 1 ne]);



%% Save dictionary

dict = struct;

dict.t2 = T2hi;
dict.b1 = b1hi;
dict.S  = Sdict_hi_norm;

save data/dict_20230629 dict


%% Now create a dictionary that is resolved over T1s

Sdict = abs(S(:,:,:,teidx));

[t2mesh,b1mesh,t1mesh] = meshgrid(T2,b1fact,T1);

%%% Up-sample this dictionary to make more values, so can just do quantized dictionary lookup

%%% upsample to 1ms in T2, 0.01 B1rel accuracy, 10ms in T1
T2hi = min(T2(:)):1:max(T2(:));
b1hi = min(b1fact(:)):0.01:max(b1fact(:));
T1hi = min(T1(:)):10:max(T1(:));

nt1hi = length(T1hi);
nt2hi = length(T2hi);
nb1hi = length(b1hi);

[t2mesh_hi,b1mesh_hi,t1mesh_hi] = meshgrid(T2hi,b1hi,T1hi);

%%% Upsample dictionary
Sdict_t1_hi = zeros([nb1hi nt2hi nt1hi ne]);
for ii=1:ne
    Sdict_t1_hi(:,:,:,ii) = interp3(t2mesh,b1mesh,t1mesh,Sdict(:,:,:,ii),t2mesh_hi,b1mesh_hi,t1mesh_hi);
end

%%% normalise the dictionary through the echoes
Snorm = vecnorm(Sdict_t1_hi,2,4); % 2-norm in 3rd direction
Sdict_t1_hi_norm = Sdict_t1_hi ./ repmat(Snorm,[1 1 1 ne]);

%% Generate dictionary - incl T1

dict = struct;

dict.t2 = T2hi;
dict.b1 = b1hi;
dict.t1 = T1hi;
dict.S  = Sdict_t1_hi_norm;

save data/dict_incT1_20230629 dict


%% Examine the error associated with fixing the T1 at 2.6s

if ~exist('dict_fixt1','var')|~exist('dict_hi','var')
    %%% load in both
    load data/dict_20230629.mat
    dict_fixt1 = dict;
    load data/dict_incT1_20230629.mat
    dict_hi=dict;

    nb1hi = length(dict_hi.b1);
    nt2hi = length(dict_hi.t2);
    nt1hi = length(dict_hi.t1);
end

%%% randomly select some samples from the full dictionary
N = 200000;
ix = fix(nb1hi*nt2hi*nt1hi*rand([1 N]));


t1_rand = zeros([N 1]);
b1_rand = zeros([N 1]);
t2m_rand = zeros([N 1]); % 'measured' t2
t2_rand = zeros([N 1]); % 'actual' t2

for ii=1:N
    [a,b,c] = ind2sub([nb1hi nt2hi nt1hi],ix(ii));

    % pull out the entry from the full dictionary
    x = squeeze(dict_hi.S(a,b,c,:));

    % match with the fixed T1 dictionary
    t2m_rand(ii) = dict_match(x,dict_hi.b1(a),dict_fixt1);

    % save actual params
    t1_rand(ii) = dict_hi.t1(c);
    t2_rand(ii) = dict_hi.t2(b);
    b1_rand(ii) = dict_hi.b1(a);
end




%% Look at statistics that are binned by B1 range...

% sort the data in increasing T1
[t1s,ix] = sort(t1_rand);
t2s = t2_rand(ix);
b1s = b1_rand(ix);
t2m = t2m_rand(ix);
t2e = 100*(t2m(:)-t2s(:))./t2s(:);

t1u = unique(t1s);

%%% which B1 ranges to look at
b1_bins = 0.2:0.05:1.5;
nbinb1 = length(b1_bins)-1;

t2errmed = zeros([length(t1u) nbinb1]);

for ii=1:length(t1u)
    ix = find(t1s==t1u(ii));
    tmp = t2e(ix);
    tmpb1 = b1s(ix);

    for jj=1:nbinb1
        ix2 = find((tmpb1>b1_bins(jj))&(tmpb1<b1_bins(jj+1)));
        tmp2 = tmp(ix2);
        t2errmed(ii,jj) = mean(tmp2);
    end
end

figfp(1)
imagesc(b1_bins(1:end-1),t1u,t2errmed,[-20 20])
cm = colorcet('d1');
colormap(cm)
grid on
hold on
[cc,ll] = contour(b1_bins(1:end-1),t1u,...
    medfilt2(t2errmed,[8 3]),[-20 -15 -10:2:-2 -1:1 2 4],'ShowText','on');
ylim([600 3800])
ll.EdgeColor = [1 1 1];
ll.LineWidth = 2;

set(gca,'FontSize',13)
title('Mean T_2 estimation error')
axis xy
xlabel('B_1^{rel}')
ylabel('T_1 (ms)')
cc=colorbar;
cc.Label.String = 'T_2 error (%)';

setpospap([400 400 600 550])
%  print -dpng -r300 t2dict_t1err.png

%% Compare the low res and upsampled dictionaries
%% Visualise the interpolated and original dictionaries

figfp(3)
nr = 3;
nc = 2;
for ii = 1:3
    subplot(nr,nc,2*(ii-1)+1)
    imagesc(T2,b1fact,Sdict_fixT1(:,:,ii),[0 7000])
    xlabel('T_2 (ms)')
    ylabel('B_1^{rel}')
    title(sprintf('Echo time #%d',ii))

    subplot(nr,nc,2*(ii-1)+2)
    imagesc(T2hi,b1hi,Sdict_hi(:,:,ii),[0 7000])
    xlabel('T_2 (ms)')
    ylabel('B_1^{rel}')
    title(sprintf('Echo time #%d',ii))
end

setpospap([400 80 550 720])