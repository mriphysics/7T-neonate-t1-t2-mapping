%% reposition subplots and switch off axes

function pp = spsjm(varargin)

%% args
gap=0.01;
top_gap=0;
left_gap=0;
right_gap=0;
bottom_gap=0;
axis_off=true;
auto=false;
for ii=1:length(varargin)
    if strcmpi(varargin{ii},'gap')
        gap = varargin{ii+1};
    end
    % ..
    if strcmpi(varargin{ii},'idx') % order of creation
        idx = varargin{ii+1};
        idx = fliplr(idx);
    end
    % gap at top
    if strcmpi(varargin{ii},'tgap')
        top_gap = varargin{ii+1};
    end
    % gap at left
    if strcmpi(varargin{ii},'lgap')
        left_gap = varargin{ii+1};
    end
    % gap at right
    if strcmpi(varargin{ii},'rgap')
        right_gap = varargin{ii+1};
    end
    % gap at bottom
    if strcmpi(varargin{ii},'bgap')
        bottom_gap = varargin{ii+1};
    end
    % do not turn off axes
    if strcmpi(varargin{ii},'axon')
        axis_off=false;
    end
    % auto => automatically order plots to have indices s.t. 1=top left etc
    if strcmpi(varargin{ii},'auto')
        auto=true;
    end
    % make all gaps the same
    if strcmpi(varargin{ii},'allgaps')
        gp = varargin{ii+1};
        bottom_gap=gp;
        right_gap=gp;
        left_gap=gp;
        top_gap=gp;
        gap=gp;
    end
end

%% get subplot handles
pp = get(gcf,'children');
% % reverse this
% pp = flipud(pp);

% 14-4-10: colorbar can interrupt this, should ignore colorbars?
cbidx=[];
for ii=1:length(pp)
    if strcmpi(get(pp(ii),'Tag'),'Colorbar')
        cbidx = [cbidx ii];
    end
    % 13-6-11: Problem with menu titles etc (?) also dump these
    if strncmp(get(pp(ii),'Tag'),'fig',3)||strncmp(get(pp(ii),'Tag'),'Fig',3)
        cbidx = [cbidx ii];
    end
    % 28-6-11: not visible entries (not sure what they are) also dumped
    if strcmpi(get(pp(ii),'Visible'),'off')
        cbidx = [cbidx ii];
    end
    % 5-10-11: Also dump legends
    if strcmpi(get(pp(ii),'Tag'),'legend')
        cbidx = [cbidx ii];
    end

    % 15-10-24: annotation
    if strcmpi(get(pp(ii),'Tag'),'scribeOverlay')
        cbidx = [cbidx ii];
    end
end
% just dump colorbar entries
pp(cbidx)=[];

%=== 19-1-11 - auto re-order
if ~auto
    % Default is reverse
    if ~exist('idx','var')
        idx = fliplr(1:length(pp));
    end
    pp(idx)=pp;
    % look at positions
    pos=zeros(length(pp),4);
    for ii=1:length(pp)
        pos(ii,:) = get(pp(ii),'position');
    end
    % #rows/cols = #unique values of the x/y coords
    nc = length(unique(pos(:,1)));
    nr = length(unique(pos(:,2)));

else
    pos=zeros(length(pp),4);
    for ii=1:length(pp)
        pos(ii,:) = get(pp(ii),'position');
    end
    posint=fix(1e2*pos); % integer representation for avoiding comparing floats
    % #rows/cols = #unique values of the x/y coords
    colpos=sort(unique(posint(:,1)),'ascend'); % make integers for 'unique'
    rowpos=sort(unique(posint(:,2)),'descend');
    nc = length(colpos);
    nr = length(rowpos);
    np = length(pp);
    % sort by position
    cidx=zeros([np 1]);
    ridx=zeros([np 1]);
    for ii=1:np
        cidx(ii) = find(posint(ii,1)==colpos);
        ridx(ii) = find(posint(ii,2)==rowpos);
    end
    % get linear index from this
    idx = nc*(ridx-1)+cidx;
    % sort
    pp(idx)=pp;
end

%% the meat of it
if numel(gap)==1
    gap = [gap gap];
end
% reposition
wc=(1-left_gap-right_gap-(nc-1)*gap(1))/nc;
wr=(1-top_gap-bottom_gap-(nr-1)*gap(2))/nr;
for ii=1:nr
    for jj=1:nc
        idx = jj + nc*(ii-1);
        if idx>length(pp)
            break
        end
        set(pp(idx),'position',[left_gap+(wc+gap(1))*(jj-1) bottom_gap+(gap(2)+wr)*(nr-ii) wc wr]);
        if axis_off
            axis(pp(idx),'off')
        end
   end
end
%%
end