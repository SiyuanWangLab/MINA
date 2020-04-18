function yy = randomizeCompartments(coeff)
% randomize compartment assignment but maintain the number of TADs in Compartments A and B

NumA = length(find(coeff>0));
NumB = length(find(coeff<=0));

Nums = 1:(NumA+NumB);
RandA = [];
rng('default');
for i = 1:NumA
    flag = 0;
    while flag == 0
        RandNum = ceil(rand*(NumA+NumB));
        if length(find(Nums==RandNum))==1
            flag = 1;
            RandA = [RandA RandNum];
            Nums = Nums(find(Nums~=RandNum));
        end
    end
end
RandB = Nums;
Compartments = [find(coeff>0); find(coeff<=0)];
CompartmentA_Random = Compartments(RandA);
CompartmentB_Random = Compartments(RandB);
yy = ones(size(coeff));
yy(CompartmentB_Random) = -1;


