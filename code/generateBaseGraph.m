function ret = generateBaseGraph()

x = [0:0.01:4];
figure(1)
plot(x, baseFun(x,1))
hold on
plot(x, baseFun(x,2), '--')
plot(x, baseFun(x,3), ':')
hold off
%set(gca, 'XTick', 1:3, 'XTickLabel', {'x_1', 'x_2', 'x_3'})
ylabel('\phi_i(x)', 'Rotation', 90, 'Fontsize', 11)


xtl = {'x_1' 'x_2' 'x_3'};
h = my_xticklabels(gca,[1 2 3],xtl, 'Fontsize', 11);
% vertical
%h = my_xticklabels([1 2 3],xtl, ...
%    'Rotation',-90, ...
%    'VerticalAlignment','middle', ...
%    'HorizontalAlignment','left');

function phi = baseFun(x, id)

phi = zeros(size(x,1), size(x,2));

low = (x >=(id-1)) & (x < id);
high = (x >= id) & (x <= (id+1));

ilow = find(low);
ihigh = find(high);

phi(ilow) = x(ilow)-id+1;
phi(ihigh) = id-x(ihigh)+1;


function ht = my_xticklabels(varargin)

% ht = my_xticklabels(Ha, xtickpos, xtickstring)
% or
% ht = my_xticklabels(xtickpos, xtickstring)

% Pekka Kumpulainen 12.2.2008

textopts = {};
if length(varargin{1})==1 && ...
        ishandle(varargin{1}) && ...
        strcmpi(get(varargin{1},'Type'),'axes');
    Ha = varargin{1};
    xtickpos = varargin{2};
    xtickstring = varargin{3};
    if nargin > 3
        textopts = varargin(4:end);
    end
else
    Ha = gca;
    Hfig = get(Ha,'Parent');
    xtickpos = varargin{1};
    xtickstring = varargin{2};
    if nargin > 2
        textopts = varargin(3:end);
    end
end

set(Ha,'XTick',xtickpos, 'XTickLabel','')
h_olds = findobj(Ha, 'Tag', 'MUXTL');
if ~isempty(h_olds)
    delete(h_olds)
end

%% Make XTickLabels 
NTick = length(xtickpos);
Ybot = min(get(gca,'YLim'));
ht = zeros(NTick,1);
for ii = 1:NTick
    ht(ii) = text('String',xtickstring{ii}, ...
        'Units','data', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center ', ...
        'Position',[xtickpos(ii) Ybot], ...
        'Tag','MUXTL');
end
if ~isempty(textopts)
    set(ht,textopts{:})
end

%% squeeze axis if needed

set(Ha,'Units','pixels')
Axpos = get(Ha,'Position');
% set(Hfig,'Units','pixels')
% Figpos = get(Hfig,'Position');

set(ht,'Units','pixels')
TickExt = zeros(NTick,4);
for ii = 1:NTick
    TickExt(ii,:) = get(ht(ii),'Extent');
end

needmove = -(Axpos(2) + min(TickExt(:,2)));

if needmove>0;
    Axpos(2) = Axpos(2)+needmove+2;
    Axpos(4) = Axpos(4)-needmove+2;
    set(Ha,'Position',Axpos);
end

set(Ha,'Units','normalized')
set(ht,'Units','normalized')