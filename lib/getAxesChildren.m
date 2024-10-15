function gg_out = getAxesChildren(h)
% function gg = getAxesChildren(h)
% gets all children of a function handle that are either axes, color bars
% or other objects. Aims to ignore objects that aren't relevant for visible
% formatting such as hidden items, context overlays etc.
% Shaihan Malik Oct 2024
gg = get(h,'Children');

% items to keep
gg_out = [];

for ii=1:length(gg),
    
    if strcmpi(gg(ii).Type,'axes')
        gg_out = [gg_out gg(ii)];
        continue
    end

    if strcmpi(gg(ii).Type,'colorbar')
        gg_out = [gg_out gg(ii)];
        continue
    end

    if strcmpi(gg(ii).Type,'legend')
        gg_out = [gg_out gg(ii)];
        continue
    end
    
end