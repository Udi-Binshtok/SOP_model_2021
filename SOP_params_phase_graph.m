function [ pattern ] = SOP_params_phase_graph( A,B,curent_project, model ,d0,d_M0 )
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

SOPs = NaN(length(valueA),length(valueB));

%lab
% CodeFolder = 'D:\Students\Udi\Udi Google Drive\University\David Spirnzak''s lab\Computer Simulation\קוד עדכני\Model for lateral inhibition of SOP 24.8.16';
%laptop
% CodeFolder = 'C:\Udi\Udi Google Drive\University\David Spirnzak''s lab\Computer Simulation\??? ?????\Model for lateral inhibition of SOP 24.8.16\2020';
% 
% cd([CodeFolder '\'])

name_of_cuurent_project = curent_project; %'lateral inhibition of SOP';
function_name = 'SOP_params_phase_graph';
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
        k = params.P*params.Q; %number of cells
        y0 = SOP_InitialConditions_rho(model.zone,params,k);
    case 'epsilon'
        params = SOP_DefaultParams_epsilon(model.zone);
        k = params.P*params.Q; %number of cells
        y0 = SOP_InitialConditions_epsilon(model.zone,params,k);
end
        

% original_valueA = params.(fieldA);
original_valueA = params.(fieldA).d; % beta.d
% original_valueA = params.(fieldB);
original_valueB = params.(fieldB).c; % Kappa.c

%%

w = waitbar(0,'calcuating');

for i = 1:length(valueA)
    
    for j = 1:length(valueB)
        %% Select which parameters to change:
        
        cl = clock; hour = num2str(cl(4)); minute = num2str(cl(5));
        fname = [FiguresFolder hour 'hr' minute 'min_' fieldA num2str(i) fieldB num2str(j)];
        
%         params.(fieldA) = valueA(i);
        params.(fieldA).d = valueA(i); % beta.d
        y0(1:k) = d0(i).*ones(k,1); % if beta.d is changed then also d0 is changed 
        y0(k+1:2*k) = d_M0(i).*ones(k,1); % if beta.d is changed then also d_M0 is changed
%         params.(fieldB) = valueA(i);
        params.(fieldB).c = valueB(j); % Kappa.c
        
        [yout,tout,params] = SOP_multicell_LI(model,params,y0);

        %% 
        
        SOP_cell = params.SOP.cell;
        SOP_nearest_neighbors = params.SOP.neares_neighbors;
        SOP_next_nearest_neighbors = params.SOP.next_nearest_neighbors;
        
        k = params.P*params.Q; %number of cells
        g.All_cells = yout(end,5*k+1:6*k);
        g.SOP_cell = g.All_cells(SOP_cell);
        g.SOP_nearest_neighbors = g.All_cells(SOP_nearest_neighbors);
        g.SOP_next_nearest_neighbors = g.All_cells(SOP_next_nearest_neighbors);
        Threshold_SOP = params.Threshold_SOP;
        fate.All_cell = min([floor(g.All_cells./Threshold_SOP) ; ones(1,k)]);
        fate.SOP_cell = fate.All_cell(SOP_cell);
        fate.SOP_nearest_neighbors = fate.All_cell(SOP_nearest_neighbors);
        fate.SOP_next_nearest_neighbors = fate.All_cell(SOP_next_nearest_neighbors);
        
        fates = [ fate.SOP_cell , fate.SOP_nearest_neighbors , fate.SOP_next_nearest_neighbors ];
        
        SOPs(i,j) = sum(fates);
        
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
        
        pattern(i,j).SOPs = sum(fates);
%         pattern(i,j).(fieldA) = params.(fieldA);
        pattern(i,j).(fieldA) = params.(fieldA); % beta.d
%         pattern(i,j).(fieldB) = params.(fieldB);
        pattern(i,j).(fieldB).c = params.(fieldB).c; % Kappa.c
        pattern(i,j).simulation.yout = yout;
        pattern(i,j).simulation.tout = tout;
        pattern(i,j).simulation.params = params;
        
        
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
    waitbar(i/length(valueA),w)
end
close(w)

% % % % ResultName = [ResultsFolder hour 'hr' minute 'min_' fieldA '_Vs_' fieldB '_params_phase_graph'];
% % % % save(ResultName,'pattern','-v7.3') % NOTE: The file is to big

% plot of SOPs pattern Vs. parameters

ff(1) = figure('NumberTitle','off','Name','stable pattern (params_phase_graph)');
surf(valueB,valueA,SOPs)
ylabel([fieldA ' [dimensionless]'],'FontSize',22)
% xlabel([fieldB ' [dimensionless]'],'FontSize',22)
xlabel([fieldB '.c [dimensionless]'],'FontSize',22)
zlabel('SOPs pattern []')

ff(2) = figure('NumberTitle','off','Name','stable pattern (params_phase_graph)');
PercentA = 100.*(valueA - original_valueA)./original_valueA;
PercentB = 100.*(valueB - original_valueB)./original_valueB;
surf(PercentB,PercentA,SOPs)
ylabel([fieldA ' [dimensionless]'],'FontSize',22)
% xlabel([fieldB ' [dimensionless]'],'FontSize',22)
xlabel([fieldB '.c [dimensionless]'],'FontSize',22)

zlabel('SOPs pattern []')

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