% 2024-09-16: voxel-wise dictionary lookup function using dot product
% matching. This version uses an input T1 as well as B1 value
% Inputs:   x - signal vector over 3 echoes
%           b1 - relative B1 value of voxel
%           t1 - t1 of voxel (ms)
%           dict = struct containing dictionary, made by dict_create.m
% Shaihan Malik, King's College London, 2024

function [t2,s0] = dict_match_inct1(x,b1,t1,dict)

%%% extract the correct B1 from the dict
[~,b1idx] = min(abs(dict.b1-b1));

%%% extract the correct B1 from the dict
[~,t1idx] = min(abs(dict.t1-t1));

%%% get T2 dictionary
D = squeeze(dict.S(b1idx,:,t1idx,:));

%%% Now dot product
dp = D * x(:);

%%% get correct T2
[s0,t2idx] = max(dp);
t2 = dict.t2(t2idx);

% % if t2 is out of range, replace with NaN
% if ismember(t2idx,[1 length(dict.t2)])
%     t2 = NaN;
% end

end
