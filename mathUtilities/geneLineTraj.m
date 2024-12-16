function traj = geneLineTraj(limits, repNum, numPerLine)
% generate traj [0,Limit], repeat repNum times.
n = length(limits);
traj = zeros(n,repNum*numPerLine+1);
    for i = 1:n
        traj(i,:) = interp1(0:repNum, limits(i)*rem(0:repNum,2), linspace(0,repNum,repNum*numPerLine+1));
    end

end