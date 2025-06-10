box off

try
    set(gca,'LineWidth',3)
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 20)
    set(xAX,'color','k')
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 20)
    set(yAX,'color','k')
    set(gca, 'TickDir', 'out')
    set(gcf,'color','w')
    set(gca,'Layer','top')
    axis normal

catch

    % Get the current figure and axis handles
    fig_handle = gcf;
    ax_handle = findobj(fig_handle, 'Type', 'axes');

    set(ax_handle,'LineWidth',3)
    xAX = get(ax_handle,'XAxis');
    set(xAX,'FontSize', 20)
    set(xAX,'color','k')
    yAX = get(ax_handle,'YAxis');
    set(yAX,'FontSize', 20)
    set(yAX,'color','k')
    set(ax_handle, 'TickDir', 'out')
    set(gcf,'color','w')
    set(ax_handle,'Layer','top')
    axis normal
end