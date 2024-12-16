function y = lowPassPadding(x, fpass, samplingRate, pad)

% add padding of the input matrix x and trim the padding.
% x: column-wise data. e.g. nx3 shape matrix 
    x_pad = [repmat(x(1,:),[pad,1]); 
             x; 
             repmat(x(end,:),[pad,1])];
    y_pad = lowpass(x_pad, fpass, samplingRate);
    y = y_pad(pad+1:end-pad, :);

end