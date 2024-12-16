function K = getTubeK(tube)

    try 
        K = tube.K;
    catch

        ns = length(tube.s);
    
        K = [tube.kb; tube.kb; tube.kt] * ones(1,ns);
        % K = diag(K(:));
    end
    K = K(:);
    
end