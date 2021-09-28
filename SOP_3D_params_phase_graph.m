function [ pattern3D, SOPs ] = SOP_3D_params_phase_graph( A,B,C,curent_project, model,d0,d_M0 )
%params_phase_graph will explore different parameters for new_changed_dimensionless
%   The function will plot out information about the stability of the
%   pattern after going through regulation process. 
%   A and B are the two parameters in test and need to be row vectors
%   containing vary values (that are already scaled as their final scailing
%   in params structure)

pp = genpath('functions');
addpath(pp)

fieldA = inputname(1);
valueA = A;
fieldB = inputname(2);
valueB = B;
fieldC = inputname(3);
valueC = C;

SOPs = NaN(length(valueA),length(valueB),length(valueC));

%lab
% CodeFolder = 'D:\Students\Udi\Udi Google Drive\University\David Spirnzak''s lab\Computer Simulation\קוד עדכני\Model for lateral inhibition of SOP 24.8.16';
%laptop
% CodeFolder = 'C:\Udi\Udi Google Drive\University\David Spirnzak''s lab\Computer Simulation\??? ?????\Model for lateral inhibition of SOP 24.8.16\2020';
% 
% cd([CodeFolder '\'])

name_of_cuurent_project = curent_project; %'lateral inhibition of SOP';
function_name = 'SOP_3D_params_phase_graph';
[ FiguresFolder, ResultsFolder] = save_to( name_of_cuurent_project, function_name , fieldA, fieldB );

% check existance of folders and create new folders if needed
if ~exist(ResultsFolder, 'dir')
    mkdir(ResultsFolder);
end
if ~exist(FiguresFolder, 'dir')
    mkdir(FiguresFolder);
end
    

% I apply params into St_multicell_LI_variable_disordered_filopodia_v1_full_equations
switch model.ratio
    case 'rho'
        params = SOP_DefaultParams_rho(model.zone);
        q = params.P*params.Q; %number of cells
        y0 = SOP_InitialConditions_rho(model.zone,params,q);
    case 'epsilon'
        params = SOP_DefaultParams_epsilon(model.zone);
        q = params.P*params.Q; %number of cells
        y0 = SOP_InitialConditions_epsilon(model.zone,params,q);
end
        

original_valueA = params.(fieldA); % rho_M_N
% original_valueA = params.(fieldA).d; % beta.d
% original_valueA = params.(fieldB);
original_valueB = params.(fieldB).c; % Kappa.c
% original_valueA = params.(fieldC);
original_valueC = params.(fieldC).d; % beta.d

%%

w = waitbar(0,'calcuating');

for i = 1:length(valueA)
    for j = 1:length(valueB)
        for k = 1:length(valueC)
        %% Select which parameters to change:
        
        cl = clock; hour = num2str(cl(4)); minute = num2str(cl(5));
        fname = [FiguresFolder hour 'hr' minute 'min_' fieldA num2str(i) fieldB num2str(j)];
        
        params.(fieldA) = valueA(i); % rho_M_N
%         params.(fieldA).d = valueA(i); % beta.d
%         params.(fieldB) = valueA(i);
        params.(fieldB).c = valueB(j); % Kappa.c
%         params.(fieldC) = valueC(k); 
        params.(fieldC).d = valueC(k); % beta.d
        y0(1:q) = d0(k).*ones(q,1); % if beta.d is changed then also d0 is changed 
        y0(q+1:2*q) = d_M0(k).*ones(q,1); % if beta.d is changed then also d_M0 is changed
        
        [yout,tout,params] = SOP_multicell_LI(model,params,y0);

        %% 
        
        SOP_cell = params.SOP.cell;
        SOP_nearest_neighbors = params.SOP.neares_neighbors;
        SOP_next_nearest_neighbors = params.SOP.next_nearest_neighbors;
        
        q = params.P*params.Q; %number of cells
        g.All_cells = yout(end,5*q+1:6*q);
        g.SOP_cell = g.All_cells(SOP_cell);
        g.SOP_nearest_neighbors = g.All_cells(SOP_nearest_neighbors);
        g.SOP_next_nearest_neighbors = g.All_cells(SOP_next_nearest_neighbors);
        Threshold_SOP = params.Threshold_SOP;
        fate.All_cell = min([floor(g.All_cells./Threshold_SOP) ; ones(1,q)]);
        fate.SOP_cell = fate.All_cell(SOP_cell);
        fate.SOP_nearest_neighbors = fate.All_cell(SOP_nearest_neighbors);
        fate.SOP_next_nearest_neighbors = fate.All_cell(SOP_next_nearest_neighbors);
        
        fates = [ fate.SOP_cell , fate.SOP_nearest_neighbors , fate.SOP_next_nearest_neighbors ];
        
        SOPs(i,j,k) = sum(fates);
        
% %         switch fates % Error: switch does not work on vectors...
% %             case zeros(1,19) 
% %                 % none of cells has been selected
% %                 SOPs(i,j) = 0;
% %                 pattern(i,j).SOPs = 0;
% %                 
% %             case [ 1 , zeros(1,18) ] 
% %                 % one SOP, at the center, has been selected
% %                 SOPs(i,j) = 1;
% %                 pattern(i,j).SOPs = 1;
% %                 
% %             case [ 0 , 1 0 1 0 1 0 , zeros(1,12) ]
% %                 % three SOPs, only among the nearest neighbors, have been selected
% %                 SOPs(i,j) = 3;
% %                 pattern(i,j).SOPs = 3;
% %                 
% %             case [ 0 , 0 1 0 1 0 1 , zeros(1,12) ]
% %                 % three SOPs, only among the nearest neighbors, have been selected
% %                 SOPs(i,j) = 3;
% %                 pattern(i,j).SOPs = 3;    
% %                 
% %             case [ 1 , zeros(1,6), 1 0 1 0 1 0 1 0 1 0 1 0 ]
% %                 % seven SOPs, one at the center and 6 next nearest neighbors, have been selected
% %                 SOPs(i,j) = 7;
% %                 pattern(i,j).SOPs = 7;
% %                 
% %             case [ 1 , zeros(1,6), 0 1 0 1 0 1 0 1 0 1 0 1 ]
% %                 % seven SOPs, one at the center and 6 next nearest neighbors, have been selected
% %                 SOPs(i,j) = 7;
% %                 pattern(i,j).SOPs = 7;
% %                 
% %             case ones(1,19)
% %                 % all of the cells have been selected
% %                 SOPs(i,j) = 19;
% %                 pattern(i,j).SOPs = 19;
% %                 
% %         end
        
        pattern3D(i,j,k).SOPs = sum(fates);
        pattern(i,j,k).(fieldA) = params.(fieldA); % rho_M_N
%         pattern3D(i,j).(fieldA).d = params.(fieldA).d; % beta.d
%         pattern(i,j).(fieldB) = params.(fieldB);
        pattern3D(i,j,k).(fieldB).c = params.(fieldB).c; % Kappa.c
%         pattern(i,j).(fieldC) = params.(fieldC);
        pattern3D(i,j,k).(fieldC).d = params.(fieldC).d; % beta.d
        pattern3D(i,j,k).simulation.yout = yout;
        pattern3D(i,j,k).simulation.tout = tout;
        pattern3D(i,j,k).simulation.params = params;
        
        
        %% Ploting

%         fa = figure('NumberTitle','off','Name',[fieldA num2str(i) fieldB num2str(j) '_a']);
%         tind = length(tout);   % final frame
%         clf;
%         for p = 1:params.P
%             for q = 1:params.Q
%                 ind = pq2ind(p,q,params.P);
% %                  mycolor = min([yout(tind,4*k+ind)/ Norm,1]); % defined the normalized color of cell ; by Repressor
%                 mycolor = min([yout(tind,k+ind)/ Norm,1]); % defined the normalized color of cell ; by Activated-Delta
%                 plotHexagon(p,q,[1,1-mycolor,1-mycolor]);
%             end
%         end
%         
%         axis image; axis off; box off;
%         
%         savefig(fa,[fname '_a'])
%         close(fa)
%         
% 
%         fb = figure('NumberTitle','off','Name',[fieldA num2str(i) fieldB num2str(j) '_b']);
%         semilogy(tout,size_diff)
%         xlabel('time [dimensionless]','FontSize',22)
%         ylabel('Absolute change in d* [ protein/time unit ]','FontSize',22)
%         
%         savefig(fb,[fname '_b'])
%         close(fb)
%         
        end
    end
    waitbar(i/length(valueA),w)
end
close(w)

% % % % ResultName = [ResultsFolder hour 'hr' minute 'min_' fieldA '_Vs_' fieldB '_params_phase_graph'];
% % % % save(ResultName,'pattern','-v7.3') % NOTE: The file is to big

% plot of SOPs pattern Vs. parameters

SOPs_default = ceil(length(valueA)/2); % Which element in SOPs is the default parameters set, assuming length(valueA) = length(valueB) = length(valueC)
PercentA = 100.*(valueA - original_valueA)./original_valueA;
PercentB = 100.*(valueB - original_valueB)./original_valueB;
PercentC = 100.*(valueC - original_valueC)./original_valueC;

%%% 2D phase diagrams
[XX,YY] = meshgrid(PercentA,PercentB);
% figure
% Kc_rho = SOPs(:,:,SOPs_default);
% surf(PercentB,PercentA,Kc_rho)
% 
% figure
% Kc_beta_T = transpose(reshape(SOPs(SOPs_default,:,:),length(valueA),length(valueA)));
% surf(PercentB,PercentC,Kc_beta_T)

figure % Kappa_c vs. rho_N^M
Kc_rho = SOPs(:,:,SOPs_default);
Kc_rho_NaN = Kc_rho;
Kc_rho_NaN(Kc_rho_NaN == 0) = NaN;
Kc_rho_NaN(Kc_rho_NaN == 1) = NaN;
% s_Kc_rho_NaN = surf(PercentB,PercentA,Kc_rho_NaN);
% s_Kc_rho_NaN.EdgeColor = 'none';
scatter(XX(:),YY(:),70,Kc_rho_NaN(:),'filled')
hold on
% Kc_rho_OnlyOne = NaN(size(Kc_rho,1),size(Kc_rho,2));
% Kc_rho_OnlyOne(find(Kc_rho == 1)) = 1;
Kc_rho_OnlyOne = find(Kc_rho == 1);
% C(:,:,1) = ones(size(Kc_rho,1),size(Kc_rho,2));
% C(:,:,2) = zeros(size(Kc_rho,1),size(Kc_rho,2));
% C(:,:,3) = zeros(size(Kc_rho,1),size(Kc_rho,2));
% s_Kc_rho_OnlyOne = surf(PercentB,PercentA,Kc_rho_OnlyOne,C);
% s_Kc_rho_OnlyOne.EdgeColor = 'none';
scatter(XX(Kc_rho_OnlyOne),YY(Kc_rho_OnlyOne),70,[1 0 0],'filled')
hold on
scatter(XX(SOPs_default,SOPs_default),YY(SOPs_default,SOPs_default),150,[0 0 0],'filled','square')
colorbar
xlabel('\DeltaK_c [%]','FontSize',14,'FontWeight','bold')
ylabel('\Delta\rho_N^M [%]','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold');

figure % Kappa_c vs. beta_d
Kc_beta_T = transpose(reshape(SOPs(SOPs_default,:,:),length(valueA),length(valueA)));
Kc_beta_T_NaN = Kc_beta_T;
Kc_beta_T_NaN(Kc_beta_T_NaN == 0) = NaN;
Kc_beta_T_NaN(Kc_beta_T_NaN == 1) = NaN;
% s_Kc_beta_T_NaN = surf(PercentB,PercentA,Kc_beta_T_NaN);
% s_Kc_beta_T_NaN.EdgeColor = 'none';
scatter(XX(:),YY(:),70,Kc_beta_T_NaN(:),'filled')
hold on
% Kc_beta_T_OnlyOne = NaN(size(Kc_rho,1),size(Kc_rho,2));
% Kc_beta_T_OnlyOne(find(Kc_beta_T == 1)) = 1;
Kc_beta_T_OnlyOne = find(Kc_beta_T == 1);
% C(:,:,1) = ones(size(Kc_beta_T,1),size(Kc_beta_T,2));
% C(:,:,2) = zeros(size(Kc_beta_T,1),size(Kc_beta_T,2));
% C(:,:,3) = zeros(size(Kc_beta_T,1),size(Kc_beta_T,2));
% s_Kc_beta_T_OnlyOne = surf(PercentB,PercentA,Kc_beta_T_OnlyOne,C);
% s_Kc_beta_T_OnlyOne.EdgeColor = 'none';
scatter(XX(Kc_beta_T_OnlyOne),YY(Kc_beta_T_OnlyOne),70,[1 0 0],'filled')
hold on
scatter(XX(SOPs_default,SOPs_default),YY(SOPs_default,SOPs_default),150,[0 0 0],'filled','square')
colorbar
xlabel('\DeltaK_c [%]','FontSize',14,'FontWeight','bold')
ylabel('\Delta\beta_d [%]','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold');


%%% 3D phase diagram
[X,Y,Z] = meshgrid(PercentA,PercentB,PercentC);
% figure
% scatter3(X(:),Y(:),Z(:),SOPs(:)*2+1,SOPs(:),'filled') % x is y and y is x why? don't know 

figure
SOPs_NaN = SOPs;
SOPs_NaN(SOPs_NaN == 0) = NaN;
SOPs_NaN(SOPs_NaN == 1) = NaN;
scatter3(X(:),Y(:),Z(:),SOPs_NaN(:)*2,SOPs_NaN(:),'filled') % x is y and y is x why? don't know
hold on
SOPs_OnlyOne = find(SOPs == 1);
scatter3(X(SOPs_OnlyOne),Y(SOPs_OnlyOne),Z(SOPs_OnlyOne),10,[1 0 0],'filled') % paint 1 SOP points in red, size 10 
hold on
scatter3(X(SOPs_default,SOPs_default,SOPs_default),Y(SOPs_default,SOPs_default,SOPs_default),Z(SOPs_default,SOPs_default,SOPs_default),50,[0 0 0],'filled','square') % paint default parameters point in black, size 20
colorbar
xlabel('\DeltaK_c [%]','FontSize',14,'FontWeight','bold')
ylabel('\Delta\rho_N^M [%]','FontSize',14,'FontWeight','bold')
zlabel('\Delta\beta_d [%]','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold');


% ff(1) = figure('NumberTitle','off','Name','stable pattern (params_phase_graph)');
% surf(valueB,valueA,SOPs)
% ylabel([fieldA ' [dimensionless]'],'FontSize',22)
% % xlabel([fieldB ' [dimensionless]'],'FontSize',22)
% xlabel([fieldB '.c [dimensionless]'],'FontSize',22)
% zlabel('SOPs pattern []')
% 
% ff(2) = figure('NumberTitle','off','Name','stable pattern (params_phase_graph)');
% PercentA = 100.*(valueA - original_valueA)./original_valueA;
% PercentB = 100.*(valueB - original_valueB)./original_valueB;
% s = surf(PercentB,PercentA,SOPs);
% ylabel([fieldA ' [dimensionless]'],'FontSize',22)
% % xlabel([fieldB ' [dimensionless]'],'FontSize',22)
% xlabel([fieldB '.c [dimensionless]'],'FontSize',22)
% 
% zlabel('SOPs pattern []')
% 
% s.EdgeColor = 'none';

hgsave(ff,[FiguresFolder hour 'hr' minute 'min_' 'final_figures_params_phase_graph'])
end



%% additional functions
function plotHexagon(p0,q0,c)

    % this function plots a hexagon centered at hex lattice coordinates p,q

    s32 = sqrt(3)/4;

    q = q0*3/4;
    p = p0*2*s32;

    if q0/2 == round(q0/2)
       p = p+s32;
    end;

    x(1) = q-.5; x(2) = q-.25; x(3) = q+.25; x(4) = q+.5; x(5) = q+.25; x(6) = q-.25;

    y(1) = p ; y(2) = p+s32; y(3) = p+s32; y(4) = p; y(5) = p-s32; y(6) = p-s32;

    c=min(c,ones(1,3));

    patch(x,y,c,'linewidth',2);
end

function ind=pq2ind(p,q, P)
    ind = p + (q-1)*P;
end