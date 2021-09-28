function [yout,tout,params] = SOP_multicell_LI(model,params,y0)

% multicell_LI simulates lateral inhibition on a hexagonal lattice. The structure params 
% contains the model parameters of the system. TOUT is a vector containing the time points 
% of the solution between 0 and Tmax. YOUT is a matrix containing the numerical solution for each variable 
% for each time point. Each row in YOUT is a vector the size of TOUT. F is a movie showing the simulation. 
%
%   model.ratio = 'rho' or 'epsilon'
%   model.zone = 'Mib1 mutual inhibition zone' or 'Neur lateral inhibition zone' or 'Gaussian-like proneural genes' 

pp = genpath('functions');
addpath(pp)

%get the default parameters, if none provided
if(nargin < 2)
    switch model.ratio
        case 'rho'
            params = SOP_DefaultParams_rho(model.zone); 
%             params = Copy_of_SOP_DefaultParams_rho(model.zone);
        case 'epsilon'
            params = SOP_DefaultParams_epsilon(model.zone);
    end
end

Tmax = params.Tmax; tspan=[0 Tmax]; %set time for simulation

% %get the connectivity matrix
params.connectivity = getconnectivityM(params.P,params.Q);
k = params.P*params.Q; %number of cells

%setting the initial conditions + noise, if none provided
if(nargin < 3)
    switch model.ratio
        case 'rho'
            y0 = SOP_InitialConditions_rho(model.zone,params,k);
%             y0 = Copy_of_SOP_InitialConditions_rho(model.zone,params,k);
        case 'epsilon'
            y0 = SOP_InitialConditions_epsilon(model.zone,params,k);
    end
end

%run simulation with lateral inhibition
[tout,yout] = ode15s(@li,tspan,y0,[],params);

% % show time traces of two cells with lateral inhibition
% 
% plot2cells(tout,yout,k)

% show lattice simulation
F=movielattice(tout,yout,k,params);

%keyboard
end 

function dy = li(t,y,params) %#ok<*INUSL>

C = params.connectivity;
k = length(C);

%%% parameters
beta_d = params.beta.d; 
beta_n = params.beta.n;
beta_E = params.beta.E;
beta_g = params.beta.g;
beta_N = params.beta.N;
beta_N_default = params.beta.N_default; % this is the Neur expression rate as in the deafult parameters 
alpha_N = params.alpha.N;
rho_M_N = params.rho_M_N;  
alpha_minus = params.alpha.minus;
Kappa_t = params.Kappa.t;
Kappa_c = params.Kappa.c;
epsilon_non = params.epsilon.non;
epsilon_M = params.epsilon.M;
T_E = params.T.E;
T_g = params.T.g;
c_s = params.c.s;
c_E = params.c.E;
c_g = params.c.g;

%%% variables
y(y < 0) = 0;   % any negative values are set to zero
d = y(1:k);         % Delta level in cells 1 to k
d_M = y(k+1:2*k);   % Mib1-activated-Delta level in cells 1 to k
d_N = y(2*k+1:3*k); % Neur-activated-Delta level in cells 1 to k
n = y(3*k+1:4*k);   % Notch level in cells 1 to k
E = y(4*k+1:5*k);   % E(spl) level in cells 1 to k
g = y(5*k+1:6*k);   % Ac and Sc level in cells 1 to k
N = y(6*k+1:7*k);   % Neur level in cells 1 to k
M = y(7*k+1:8*k);   % Mib1 level in cells 1 to k

nneighbor = C*n;    % Notch level in the neighboring cells
dneighbor = C*d;    % Delta level in the neighboring cells
d_Mneighbor = C*d_M;% Mib1-activated-Delta level in the neighboring cells 
d_Nneighbor = C*d_N;% Neur-activated-Delta level in the neighboring cells

%%% differential equations (dimensionless):
dd = beta_d + alpha_minus.*(d_M + d_N) - (alpha_N.*beta_N_default.*rho_M_N + alpha_N.*N).*d - (1 + Kappa_t.*epsilon_non.*nneighbor + Kappa_c.*n).*d; % Delta level differential equation
dd_M = alpha_N.*beta_N_default.*rho_M_N.*d - (alpha_minus + 1 + Kappa_t.*epsilon_M.*nneighbor + Kappa_c.*n).*d_M; % Mib1-activated-Delta level differential equation
dd_N = alpha_N.*N.*d - (alpha_minus + 1 + Kappa_t.*nneighbor + Kappa_c.*n).*d_N; % Neur-activated-Delta level differential equation
dn = beta_n - n - Kappa_t.*(epsilon_non.*dneighbor + epsilon_M.*d_Mneighbor + d_Nneighbor).*n - Kappa_c.*(d + d_M + d_N).*n; % Notch level differential equation
dE = beta_E.*(((Kappa_t.*(epsilon_non.*dneighbor + epsilon_M.*d_Mneighbor + d_Nneighbor).*n).^c_s)./(1+(Kappa_t.*(epsilon_non.*dneighbor + epsilon_M.*d_Mneighbor + d_Nneighbor).*n).^c_s)) - E; % E(spl) level differential equation
dg = beta_g./(1+(E./T_E).^c_E) - g; % Ac and Sc levels differential equation
dN = beta_N.*(((g./T_g).^c_g)./(1 + (g./T_g).^c_g)) - N; % Neur level differential equation
dM = zeros(k,1); % Mib1 level differential equation

dy = [dd;dd_M;dd_N;dn;dE;dg;dN;dM];

end

function C=getconnectivityM(P,Q)

k=P*Q; %number of cells
C=zeros(k,k); %This is the connectivity matrix
w=1/6; % Weight for interactions

% calculating the connectivity matrix
for s=1:k
    kneighbour=findneighbourhex(s,P,Q); %finds the neighbors of cell s
    for r=1:6
        C(s,kneighbour(r))=w;
    end
end

end 

function plot2cells(tout,yout,k)

figure(21)
clf
for i=1:2
    subplot(1,2,i)
    plot(tout,yout(:,i),'-r','linewidth',2)   %plot D levels 
    hold on
    plot(tout,yout(:,k+i),'-b','linewidth',2) %plot R levels
    title(['cell #',num2str(i)])
    xlabel('normalized time'); ylabel('normalized concentration')
    legend('D','R')
end

end

function out = findneighbourhex(ind,P,Q)
[p,q] = ind2pq(ind,P);

%above and below:
out(1) = pq2ind(mod(p,P)+1,q,P);
out(2) = pq2ind(mod(p-2,P)+1,q,P);

%left side:
qleft = mod(q-2,Q)+1;
qright = mod(q,Q)+1;

if q/2~=round(q/2)
    pup = p;
    pdown = mod(p-2,P)+1;
else 
    pup = mod(p,P)+1;
    pdown = p;
end;
out(3) = pq2ind(pup,qleft,P);
out(4) = pq2ind(pdown,qleft,P);
out(5) = pq2ind(pup,qright,P);
out(6) = pq2ind(pdown,qright,P);

end

function ind=pq2ind(p,q, P)
ind = p + (q-1)*P;
end

function [p,q]=ind2pq(ind, P)
q = 1+floor((ind-1)/P);
p = ind - (q-1)*P;
end

function plotHexagon(p0,q0,c)

% this function plots a hexagon centered at hex lattice coordinates p,q

s32 = sqrt(3)/4;

q = q0*3/4;
p = p0*2*s32;

if q0/2 == round(q0/2),
   p = p+s32;
end;

x(1) = q-.5; x(2) = q-.25; x(3) = q+.25; x(4) = q+.5; x(5) = q+.25; x(6) = q-.25;

y(1) = p ; y(2) = p+s32; y(3) = p+s32; y(4) = p; y(5) = p-s32; y(6) = p-s32;

c=min(c,ones(1,3));

patch(x,y,c,'linewidth',2);

end

function F=movielattice(tout,yout,k,params)

Fig = figure(23);
P = params.P;
Q = params.Q;
Threshold_Neur = params.T.g;
Threshold_SOP = params.Threshold_SOP;
% frameind=0;
FrameIndex = 1;
NumOfFrames = 1; %10; % number of frames for each time step;
MovieName = 'movie_SOP_model';
    
for tind = unique(round(linspace(1,length(tout),10))) % show 10 frames
% for tind = 1 % show only initial Ac/Sc pattern 
    clf;
    for i = 1:P
        for j = 1:Q
            ind = pq2ind(i,j,P);
            % define cell color by Ac and Sc level reletive to threshold
%             mycolor = min([yout(tind,5*k+ind)./Threshold_SOP,1]); % present Ac and Sc 
%             plotHexagon(i,j,[1-mycolor,1-mycolor,1]); % present Ac and Sc expression in blue 
%             mycolor = min([floor(yout(tind,5*k+ind)./Threshold_SOP),1]); % present SOPs 
%             plotHexagon(i,j,[1,1-mycolor,1-mycolor]); % present SOPs in red 
            mycolor_AcSc = min([yout(tind,5*k+ind)./Threshold_Neur,1]); % present Ac and Sc + SOPs
            mycolor_SOP = min([floor(yout(tind,5*k+ind)./Threshold_SOP),1]);
            plotHexagon(i,j,[1-mycolor_AcSc*(1-mycolor_SOP),1-mycolor_AcSc,1*(1-mycolor_SOP)]); % present Ac and Sc + SOPs  
        end
    end
    axis image; axis off; box off;
    
%     frameind=frameind+1;
%     F(frameind) = getframe; %generates a movie variable
    F = 1;
    
    %%% creating a movie
    get = getframe(Fig);
    for  m = FrameIndex:FrameIndex+NumOfFrames-1
        Movie(m) = get; %generates a movie variable
    end
    FrameIndex = FrameIndex + NumOfFrames;
%     close(Fig)
    
end

% v = VideoWriter(MovieName,'Uncompressed AVI');
% open(v)
% writeVideo(v,Movie); % save movie in avi format
% close(v)

% % if ishandle(24)
% %     close(24)
% % end
% % figure(24)
% %     timeEnd = max(find(tout<0.01*max(tout)));
% %     d_star_total = sum(yout(1:timeEnd,k+1:2*k),2)/max(sum(yout(1:timeEnd,k+1:2*k),2));
% %     d_total = sum(yout(1:timeEnd,2*k+1:3*k),2)/max(sum(yout(1:timeEnd,2*k+1:3*k),2));
% %     plot(tout(1:timeEnd),d_star_total,'r','LineWidth',6)
% %     hold on
% %     plot(tout(1:timeEnd),d_total,'-g','LineWidth',6);
% %     legend({'d^*','d'},'FontSize',30,'FontWeight','bold');
% %     xlabel('Time [a.u]','FontSize',28,'FontWeight','bold');
% %     ylabel('Normalized summation of level [dimensionless]','FontSize',28,'FontWeight','bold');
%     axes('FontSize',24,'FontWeight','bold');
% movie2avi(F,'movielattice'); % save movie in avi format
end
