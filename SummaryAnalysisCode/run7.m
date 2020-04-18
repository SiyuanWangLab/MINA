clear all
close all
% Plot the polarization of sample individual chromosomes for all cell types

load('run6workspace.mat');
mkdir figures_singleChromosomes
for ii = 1:max(CellTypeID_new)
    ChrInGroup = Chr(find([Chr.CellType]==ii));
    coeff = AllCoeff{ii};
    for i = 1:length(ChrInGroup)
        if length(find(ChrInGroup(i).r == 0))>5
            continue
        else
            try
                ChrInGroup(i).PolarizationIndex = calculatePolarization(ChrInGroup(i),coeff);
            catch
                continue
            end
        end
        if ChrInGroup(i).PolarizationIndex > 0.8
            % plot the chromosome
            ChrChozen = ChrInGroup(i);
            CompartmentA_ori = find(coeff>0);
            CompartmentB_ori = find(coeff<0);
            
            % below are codes adopted from the 2016 paper
            CompartmentA = CompartmentA_ori(find(ChrChozen.r(CompartmentA_ori)==1));
            CompartmentB = CompartmentB_ori(find(ChrChozen.r(CompartmentB_ori)==1));
            N_B = length(CompartmentB);
            N_A = length(CompartmentA);
            ChrChozen.PosB = [ChrChozen.x(CompartmentB) ChrChozen.y(CompartmentB) ChrChozen.z(CompartmentB)];
            ChrChozen.PosA = [ChrChozen.x(CompartmentA) ChrChozen.y(CompartmentA) ChrChozen.z(CompartmentA)];
            MeanB = mean(ChrChozen.PosB, 1);
            MeanA = mean(ChrChozen.PosA, 1);
            Vector = MeanB-MeanA;
            Phi = atan(Vector(2)/Vector(1));
            if Vector(1)<0
                Phi = Phi+pi;
            end
            Theta = acos(Vector(3)/(sum(Vector.^2))^0.5);
            Trasnform1 = [cos(-Phi) sin(-Phi) 0; -sin(-Phi) cos(-Phi) 0; 0 0 1];
            Transform2 = [cos(-Theta) 0 -sin(-Theta); 0 1 0; sin(-Theta) 0 cos(-Theta)];

            MEAN = repmat((MeanA+MeanB)/2, N_B,1);
            ChrChozen.PosNewB = (ChrChozen.PosB-MEAN)*Trasnform1*Transform2;
            MEAN = repmat((MeanA+MeanB)/2, N_A,1);
            ChrChozen.PosNewA = (ChrChozen.PosA-MEAN)*Trasnform1*Transform2;
            
            Pos = [ChrChozen.x ChrChozen.y ChrChozen.z];
            MEAN = repmat((MeanA+MeanB)/2, size(Pos,1),1);
            ChrChozen.PosNew = (Pos-MEAN)*Trasnform1*Transform2;
            

            
            figure(6)

            plot(ChrChozen.PosNewB(:,1), ChrChozen.PosNewB(:,3),'.b','MarkerSize', 30);
            hold on
            plot(ChrChozen.PosNewA(:,1), ChrChozen.PosNewA(:,3),'.r','MarkerSize', 30);
            axis equal

            set(gca, 'LineWidth', 1, 'FontSize', 16);
            hold off
            savefig(['figures_singleChromosomes\AB dots for Group ' num2str(ii) '.fig'])
            
            figure(5)

            plot(ChrChozen.PosNewB(:,1), ChrChozen.PosNewB(:,3),'.b','MarkerSize', 30);
            hold on
            plot(ChrChozen.PosNewA(:,1), ChrChozen.PosNewA(:,3),'.r','MarkerSize', 30);
            axis equal

            set(gca, 'LineWidth', 1, 'FontSize', 16);
            
            [K_A, V_A] = convhull(ChrChozen.PosNewA);
            [K_B, V_B] = convhull(ChrChozen.PosNewB);
            
            for n = 1:size(K_A,1)
                P1 = K_A(n,1);
                P2 = K_A(n,2);
                P3 = K_A(n,3);
                x = [ChrChozen.PosNewA(P1,1) ChrChozen.PosNewA(P2,1) ChrChozen.PosNewA(P3,1) ChrChozen.PosNewA(P1,1)];
                z = [ChrChozen.PosNewA(P1,3) ChrChozen.PosNewA(P2,3) ChrChozen.PosNewA(P3,3) ChrChozen.PosNewA(P1,3)];
                plot(x, z, '-r', 'LineWidth', 2);
            end
            for n = 1:size(K_B,1)
                P1 = K_B(n,1);
                P2 = K_B(n,2);
                P3 = K_B(n,3);
                x = [ChrChozen.PosNewB(P1,1) ChrChozen.PosNewB(P2,1) ChrChozen.PosNewB(P3,1) ChrChozen.PosNewB(P1,1)];
                z = [ChrChozen.PosNewB(P1,3) ChrChozen.PosNewB(P2,3) ChrChozen.PosNewB(P3,3) ChrChozen.PosNewB(P1,3)];
                plot(x, z, '-b','LineWidth', 2);
            end
            hold off
            savefig(['figures_singleChromosomes\AB dots with edges for Group ' num2str(ii) '.fig'])
            
            figure(7)
            [K_A, V_A] = convhull(ChrChozen.PosNewA);
            [K_B, V_B] = convhull(ChrChozen.PosNewB);
            
            DT1 = DelaunayTri(ChrChozen.PosNewA);  % Create the tetrahedral mesh
            hullFacets1 = convexHull(DT1);       % Find the facets of the convex hull
            DT2 = DelaunayTri(ChrChozen.PosNewB);  % Create the tetrahedral mesh
            hullFacets2 = convexHull(DT2);       % Find the facets of the convex hull
            
            scatter3(ChrChozen.PosNewA(:,1), ChrChozen.PosNewA(:,2), ChrChozen.PosNewA(:,3), 'r.', ...
                'SizeData', 300); hold on;
            trisurf(hullFacets1,DT1.X(:,1),DT1.X(:,2),DT1.X(:,3),'FaceColor','r', ...
                'FaceAlpha', 0.25, ...
                'EdgeColor', 'r'); hold on;
            scatter3(ChrChozen.PosNewB(:,1), ChrChozen.PosNewB(:,2), ChrChozen.PosNewB(:,3), 'b.', ...
                'SizeData', 300); hold on;
            trisurf(hullFacets2,DT2.X(:,1),DT2.X(:,2),DT2.X(:,3),'FaceColor','b', ...
                'FaceAlpha', 0.25, ...
                'EdgeColor', 'b')
            
            view([0 0])
            PlotProp
            grid off
            axis equal

            hold off
            savefig(['figures_singleChromosomes\Conv hull for Group ' num2str(ii) '.fig'])
            break
        end
    end
end
    
