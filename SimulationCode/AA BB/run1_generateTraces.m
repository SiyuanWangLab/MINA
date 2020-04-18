clear all
close all

load('G:\Miao\Summary\AllCoeff.mat')
Scores = AllCoeff{3};

KA = 50; % A interaction energy 
KB = 50; % B interaction energy
KAS = 0; % A-surface interaction energy
KBS = 0; % B-surface interaction energy
MaxDis = 4; % maximum distance between adjacent nodes.
StepsOfSimulation = 60000;
PolymerLength = length(Scores);
ANodes = find(Scores>0);
BNodes = find(Scores<0);
figure(1)
bar(Scores);

mkdir Figs
mkdir Traces
for ii = 1:100 % 100 traces to generate
    %% set initial polymer conformation
    X = 0;
    Y = 0;
    Z = 0;
    for i = 2:PolymerLength
        flag = 0;
        while flag == 0
            Xnew = X(end)+ceil(rand*3)-2;
            Ynew = Y(end)+ceil(rand*3)-2;
            Znew = Z(end)+ceil(rand*3)-2;
            if min(((X-Xnew).^2+(Y-Ynew).^2+(Z-Znew).^2).^0.5)>0
                X = [X; Xnew];
                Y = [Y; Ynew];
                Z = [Z; Znew];
                flag = 1;
            end
        end
    end

    %%
    % simultation restriction: in each step, one node on the polymer can move
    % +1 or -1 along one of the three dimentions, to a free position not occupied 
    % by another node on the polymer. Energies are calcualted before and after the move
    % and the energy differences determine the probability of the move based on
    % boltzmann distribution. The maximum allowed distance between two adjacent nodes
    % on the polymer is MaxDis. 

    n = 1; % simulation step count
    Energy(n) = calculateEnergy(X,Y,Z,Scores,ANodes,BNodes,KA,KB,KAS,KBS);
    for n = 2:StepsOfSimulation+2
        % find a movable node and an available position nearby
        flag = 0;
        while flag == 0
            N = ceil(rand*PolymerLength); % node that we attempt to move
            D = ceil(rand*3); % dimension of the movement
            if D == 1
                NewPos = [X(N)+(round(rand)-0.5)*2,Y(N),Z(N)];
            elseif D == 2
                NewPos = [X(N),Y(N)+(round(rand)-0.5)*2,Z(N)];
            else
                NewPos = [X(N),Y(N),Z(N)+(round(rand)-0.5)*2];
            end
            % determine if the NewPos is available
            if min(((X-NewPos(1)).^2+(Y-NewPos(2)).^2+(Z-NewPos(3)).^2).^0.5)>0 
                % determin if the adjacent nodes are not exceeding the MaxDis
                if N == 1
                    if sum((NewPos-[X(2) Y(2) Z(2)]).^2)>MaxDis^2
                        continue
                    end
                elseif N == PolymerLength
                    if sum((NewPos-[X(PolymerLength-1) Y(PolymerLength-1) Z(PolymerLength-1)]).^2)>MaxDis^2
                        continue
                    end
                else
                    if sum((NewPos-[X(N-1) Y(N-1) Z(N-1)]).^2)>MaxDis^2 || ...
                            sum((NewPos-[X(N+1) Y(N+1) Z(N+1)]).^2)>MaxDis^2
                        continue
                    end
                end
                % attempt to move based on energy difference. If successfully
                % moved, change flag value to 1
                Xnew = X;
                Ynew = Y;
                Znew = Z;
                Xnew(N) = NewPos(1);
                Ynew(N) = NewPos(2);
                Znew(N) = NewPos(3);
                EnergyNew = calculateEnergy(Xnew,Ynew,Znew,Scores,ANodes,BNodes,KA,KB,KAS,KBS);
                if rand<exp(Energy(end)-EnergyNew)
                    flag = 1;
                    % make transition
                    X = Xnew;
                    Y = Ynew;
                    Z = Znew;
                    Energy(n) = EnergyNew;
                end
            end
        end

    end
    %%
    figure(2)
    plot(Energy,'*b')
    savefig(['Figs/Energy' num2str(ii) '.fig'])

    figure(3)
    plot3(X,Y,Z,'-k')
    hold on
    plot3(X(ANodes),Y(ANodes),Z(ANodes),'*r')
    plot3(X(BNodes),Y(BNodes),Z(BNodes),'*b')
    axis equal
    hold off
    savefig(['Figs/Trace' num2str(ii) '.fig'])
    
    save(['Traces/Trace' num2str(ii) '.mat'], 'Energy', 'X', 'Y', 'Z');
    
    clear Energy
    clear X
    clear Y
    clear Z
end