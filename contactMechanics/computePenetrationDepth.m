function depth = computePenetrationDepth(p, plane)
% p: tube centerline

    depth = -plane.n' * (p - plane.p0);

end