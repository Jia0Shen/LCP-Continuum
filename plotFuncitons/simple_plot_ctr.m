function names = simple_plot_ctr(g_arrays, r, color)
if nargin == 1
    r = 1;
    name1 = plot3(g_arrays{1}(:,13),g_arrays{1}(:,14),g_arrays{1}(:,15), 'Color', [1 1 1]*0.8, 'LineWidth', 5*r);
    name2 = plot3(g_arrays{2}(:,13),g_arrays{2}(:,14),g_arrays{2}(:,15), 'Color', [1 1 1]*0.6, 'LineWidth', 3*r);
    name3 = plot3(g_arrays{3}(:,13),g_arrays{3}(:,14),g_arrays{3}(:,15), 'Color', [1 1 1]*0.3, 'LineWidth', r);
elseif nargin == 2
    name1 = plot3(g_arrays{1}(:,13),g_arrays{1}(:,14),g_arrays{1}(:,15), 'Color', [1 1 1]*0.8, 'LineWidth', 5*r);
    name2 = plot3(g_arrays{2}(:,13),g_arrays{2}(:,14),g_arrays{2}(:,15), 'Color', [1 1 1]*0.6, 'LineWidth', 3*r);
    name3 = plot3(g_arrays{3}(:,13),g_arrays{3}(:,14),g_arrays{3}(:,15), 'Color', [1 1 1]*0.3, 'LineWidth', r);
elseif nargin == 3
    name1 = plot3(g_arrays{1}(:,13),g_arrays{1}(:,14),g_arrays{1}(:,15), 'Color', color, 'LineWidth', 5*r);
    name2 = plot3(g_arrays{2}(:,13),g_arrays{2}(:,14),g_arrays{2}(:,15), 'Color', color, 'LineWidth', 3*r);
    name3 = plot3(g_arrays{3}(:,13),g_arrays{3}(:,14),g_arrays{3}(:,15), 'Color', color, 'LineWidth', r);
end

names = [name1, name2, name3];
end