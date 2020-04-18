function yy = calculatePolarization(ChrChozen,coeff)
CompartmentA_ori = find(coeff>0);
CompartmentB_ori = find(coeff<=0);
CompartmentA = CompartmentA_ori(find(ChrChozen.r(CompartmentA_ori)==1));
CompartmentB = CompartmentB_ori(find(ChrChozen.r(CompartmentB_ori)==1));
N_B = length(CompartmentB);
N_A = length(CompartmentA);

ChrChozen.PosB = [ChrChozen.x(CompartmentB) ChrChozen.y(CompartmentB) ChrChozen.z(CompartmentB)];
ChrChozen.PosA = [ChrChozen.x(CompartmentA) ChrChozen.y(CompartmentA) ChrChozen.z(CompartmentA)];

% Here introduce the convex hull analysis
[K_A, V_A] = convhull(ChrChozen.PosA);
[K_B, V_B] = convhull(ChrChozen.PosB);
    
Interval = 0.01; %um
PA = ChrChozen.PosA;
for i = 1:size(K_A,1)
    P1 = K_A(i,1);
    P2 = K_A(i,2);
    P3 = K_A(i,3);
    x = ChrChozen.PosA(P1,:);
    y = ChrChozen.PosA(P2,:);
    z = ChrChozen.PosA(P3,:);
    D = (sum((x-y).^2))^0.5;
    J = floor(D/Interval);
    for j = 1:J
        NewPoint = x + (y-x)*j/(D/Interval);
        PA = [PA; NewPoint];
    end
    D = (sum((y-z).^2))^0.5;
    J = floor(D/Interval);
    for j = 1:J
        NewPoint = y + (z-y)*j/(D/Interval);
        PA = [PA; NewPoint];
    end
    D = (sum((z-x).^2))^0.5;
    J = floor(D/Interval);
    for j = 1:J
        NewPoint = z + (x-z)*j/(D/Interval);
        PA = [PA; NewPoint];
    end
end
PB = ChrChozen.PosB;
for i = 1:size(K_B,1)
    P1 = K_B(i,1);
    P2 = K_B(i,2);
    P3 = K_B(i,3);
    x = ChrChozen.PosB(P1,:);
    y = ChrChozen.PosB(P2,:);
    z = ChrChozen.PosB(P3,:);
    D = (sum((x-y).^2))^0.5;
    J = floor(D/Interval);
    for j = 1:J
        NewPoint = x + (y-x)*j/(D/Interval);
        PB = [PB; NewPoint];
    end
    D = (sum((y-z).^2))^0.5;
    J = floor(D/Interval);
    for j = 1:J
        NewPoint = y + (z-y)*j/(D/Interval);
        PB = [PB; NewPoint];
    end
    D = (sum((z-x).^2))^0.5;
    J = floor(D/Interval);
    for j = 1:J
        NewPoint = z + (x-z)*j/(D/Interval);
        PB = [PB; NewPoint];
    end
end
        
AinB = [];
for i = 1:size(PA,1)
    P = [ChrChozen.PosB; PA(i,:)];
    [K_P, V_P] = convhull(P);
    if V_P <= V_B
        AinB = [AinB; PA(i,:)];
    end
end
BinA = [];
for i = 1:size(PB,1)
    P = [ChrChozen.PosA; PB(i,:)];
    [K_P, V_P] = convhull(P);
    if V_P <= V_A
        BinA = [BinA; PB(i,:)];
    end
end
if size(AinB,1)+size(BinA,1) == 0
    yy = 1;
elseif size(AinB,1)+size(BinA,1) < 4
    yy = 1;
else
    PosShared = [AinB; BinA];
    [K_Shared, V_Shared] = convhull(PosShared);
    yy = ((1-V_Shared/V_A)*(1-V_Shared/V_B))^0.5;
end