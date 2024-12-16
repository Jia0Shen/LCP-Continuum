function vid = makeVid2D(axis_handle, path, sampleSkip)

if nargin == 2
    sampleSkip = 1;
end
% make a list of plot handles to videos
h_list = findobj(axis_handle, 'Type', 'line');
h_circle = copyobj(h_list, axis_handle);
h_copy = copy(h_list);

dim_list = cellfun(@length,{h_list.XData});
max_dim = max(dim_list);
frameSample = 1:round(sampleSkip):max_dim;
n_frames = length(frameSample);
n = length(h_list);

totalTime = 10;  % 10s
vid = VideoWriter(path,'MPEG-4');
vid.FrameRate = round(n_frames/totalTime);
open(vid)

for it = frameSample

    % change the values
    for i = 1:n
        hi = h_list(i);
        hi_circle = h_circle(i);
        hi_copy = h_copy(i);

        max_it = length(hi_copy.XData);
        idx = min(it, max_it);

        set(hi, 'XData', hi_copy.XData(1:idx), ...
                'YData', hi_copy.YData(1:idx));

        % add circle at the end.
        set(hi_circle, 'XData', hi_copy.XData(idx), ...
                       'YData', hi_copy.YData(idx),...
                       'Marker', 'o');
        if it == frameSample(end)
            set(hi_circle, 'Marker', 'none')
        end

    end
  

    % save video
    drawnow
    writeVideo(vid,getframe(gcf));
    
end

close(vid)

end