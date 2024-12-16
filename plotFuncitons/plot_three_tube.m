function plot_three_tube(g_arrays, color)
    if nargin ==1
       plotSection(g_arrays{1}, 0.0015, [1 1 1]*0.4);
       plotSection(g_arrays{2}, 0.001, [1 1 1]*0.55);
       plotSection(g_arrays{3}, 0.0005, [1 1 1]*0.7);
    elseif nargin == 2
       plotSection(g_arrays{1}, 0.0015, color);
       plotSection(g_arrays{2}, 0.001, color);
       plotSection(g_arrays{3}, 0.0005, color);
    end
end