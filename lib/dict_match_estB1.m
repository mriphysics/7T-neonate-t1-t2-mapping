% 2024-09-16: voxel-wise dictionary lookup function using dot product
% matching
% Inputs:   x - signal vector over 3 echoes
%           dict = struct containing dictionary, made by dict_create.m
% This version estimates t2 and b1, as in the original Ben Eliezer method
% Shaihan Malik, King's College London, 2024

function [t2,b1,s0] = dict_match_estB1(x,dict)

%%% extract the correct B1 from the dict
% [~,b1idx] = min(abs(dict.b1-b1));

%%% get T2 dictionary
[n1, n2, n3] = size(dict.S);

D = reshape(dict.S,[n1*n2 n3]);

%%% Now dot product
dp = D * x(:);

%%% get index of max
[s0,maxidx] = max(dp);
[b1idx t2idx] = ind2sub([n1 n2],maxidx);

%%% get correct T2
t2 = dict.t2(t2idx);
b1 = dict.b1(b1idx);

end
