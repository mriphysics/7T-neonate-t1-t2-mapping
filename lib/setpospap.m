%%% 4-10-11: 
% function [] = setpospap(pos)
% Sets position to vector pos and paperpositionmode to auto
% 2023-09-04: 
function [] = setpospap(pos,color)

if ~exist('pos','var')
    pos = get(gcf,'position');
end
if nargin<2
    color = [1 1 1];
end

set(gcf,'paperPositionMode','auto',...
    'position',pos,'color',color,'InvertHardcopy','off');

end