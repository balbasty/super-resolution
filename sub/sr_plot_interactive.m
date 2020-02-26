function sr_plot_interactive(action,arg1,arg2,arg3)
% Visual output for function minimisation
%
% FORMAT sr_plot_interactive('Init', Title, LabelY, LabelX)
% FORMAT sr_plot_interactive('Set', x, y)
% FORMAT sr_plot_interactive('Clear')
persistent min1dplot

if ~nargin, action = 'Init'; end

% Find the Interactive window and exit if not
%--------------------------------------------------------------------------
fg = spm_figure('FindWin','Interactive');
if isempty(fg), return; end

%-Initialize
%--------------------------------------------------------------------------
if strcmpi(action,'init')
    if nargin<4, arg3 = 'Function';          end
    if nargin<3, arg2 = 'Value';             end
    if nargin<2, arg1 = 'Line minimisation'; end
    
    min1dplot = struct('pointer',get(fg,'Pointer'),...
                       'name',   get(fg,'Name'),...
                       'ax',     []);
    sr_plot_interactive('Clear');
    set(fg,'Pointer','Watch');
    min1dplot.ax = axes('Position', [0.15 0.1 0.8 0.75],...
                        'Box',      'on',...
                        'Parent',   fg);
    lab = get(min1dplot.ax,'Xlabel');
    set(lab,'string',arg3,'FontSize',10);
    lab = get(min1dplot.ax,'Ylabel');
    set(lab,'string',arg2,'FontSize',10);
    lab = get(min1dplot.ax,'Title');
    set(lab,'string',arg1);
    line('Xdata',[], 'Ydata',[],...
        'LineWidth',2,'Tag','LinMinPlot',...
        'LineStyle','-','Marker','o',...
        'Parent',min1dplot.ax);
    drawnow;
    
%-Reset
%--------------------------------------------------------------------------
elseif strcmpi(action,'set')
    br = findobj(fg,'Tag','LinMinPlot');
    if ~isempty(br)
        [xd,indx] = sort([get(br,'Xdata') arg1]);
        yd = [get(br,'Ydata') arg2];
        yd = yd(indx);
        set(br,'Ydata',yd,'Xdata',xd);
        drawnow;
    end
    
%-Clear
%--------------------------------------------------------------------------
elseif strcmpi(action,'clear')
    fg = spm_figure('FindWin','Interactive');
    if isstruct(min1dplot)
        if ishandle(min1dplot.ax), delete(min1dplot.ax); end
        set(fg,'Pointer',min1dplot.pointer);
        set(fg,'Name',min1dplot.name);
    end
    spm_figure('Clear',fg);
    drawnow;
end