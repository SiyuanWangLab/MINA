function yy = calculateSurfaceIndex(Chr)
Ind = find(Chr.r == 1);
Pos = [Chr.x(Ind) Chr.y(Ind) Chr.z(Ind)];
[K, V] = convhull(Pos);
for i = 1:length(Chr.r)
    if Chr.r(i) == 0
        yy(i) = NaN;
    else
        IndNew = Ind(find(Ind~=i));
        PosNew = [Chr.x(IndNew) Chr.y(IndNew) Chr.z(IndNew)];
        [KNew, VNew] = convhull(PosNew);
        if VNew<V
            yy(i) = 1;
        else
            yy(i) = 0;
        end
    end
end