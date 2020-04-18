function yy = calculateEnergy(X,Y,Z,Scores,ANodes,BNodes,KA,KB,KAS,KBS)
% calculate the energy of the polymer

% AA interaction energy
Pos = [X Y Z];
APos = Pos(ANodes,:);
AScores = Scores(ANodes);
D = pdist(APos);
Dmatrix = squareform(D);
Dmatrix(find(Dmatrix>=2))=0;
Dmatrix(find(Dmatrix<2 & Dmatrix>0))=1;
yy = -KA*(AScores'*Dmatrix*AScores)/2;

% BB interaction energy
BPos = Pos(BNodes,:);
BScores = Scores(BNodes);
D = pdist(BPos);
Dmatrix = squareform(D);
Dmatrix(find(Dmatrix>=2))=0;
Dmatrix(find(Dmatrix<2 & Dmatrix>0))=1;
yy = yy -KB*(BScores'*Dmatrix*BScores)/2;

% A-surface interaction energy
if KAS~=0
    % calculate the 3D convex hull
    K = convhull(Pos);
    for i = 1:length(ANodes)
        if length(find(K == ANodes(i)))>0
            yy = yy -KAS*AScores(i);
        end
    end
end

% B-surface interaction energy
if KBS~=0
    % calculate the 3D convex hull
    K = convhull(Pos);
    for i = 1:length(BNodes)
        if length(find(K == BNodes(i)))>0
            yy = yy -KBS*abs(BScores(i));
        end
    end
end








