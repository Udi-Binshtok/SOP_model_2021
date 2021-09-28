function params = SOP_DefaultParams_rho(model)
%DEFUALT_PARAMS Summary of this function goes here
%   model = 'Mib1 mutual inhibition zone' or 'Neur lateral inhibition zone' or 'Gaussian-like proneural genes'  

params.Tmax = 100;
params.P = 8;
params.Q = 8;

params.beta.d = 1*1.62; % until2020Dec23 1.62; % from equation on beta_d 1.63 % old version rho=alpha/alpha: 1*1.54; % 1.54; 
params.beta.n = 1*1.52; % until2020Dec23 1.52; % from equation on beta_d 1.49 % old version rho=alpha/alpha: 1*1.47; % 1.47;
params.beta.E = 1*1.62; % until2020Dec23 1.62;% old version rho=alpha/alpha: 1*1.54; % 1.54;
% params.beta.g = model dependent
MutualInhibition_beta_g = 0.0046;
LateralInhibition_beta_g_factor = 2*100; % 2*100; % 1.5*100;
params.beta.N = 1*1.62; % 1.62; % old version rho=alpha/alpha: 1*1.54; %4*1.54; % 1*1.54;
params.beta.N_default = 1.62; % this term is to keep the ratio rho constant while beta.N changes 
params.alpha.N = 0.8;
params.rho_M_N = 1*0.2; % 0.2
params.alpha.minus = 0.4;
params.Kappa.t = 0.8; % 0.8;
params.Kappa.c = 1*0.4; % 0.4;
params.epsilon.non = 0;
params.epsilon.M = 1;
params.T.E = 0.0051; % old version rho=alpha/alpha: 0.0019; % 0.0030;
params.T.g = 0.051; % old version rho=alpha/alpha: 0.019; % 0.014;
params.c.s = 2;
params.c.E = 3;
params.c.g = 5;

params.Threshold_SOP = 1.1*params.T.g; % 4*params.T.g;

switch model
    case 'Mib1 mutual inhibition zone' 
%         params.beta.g = 0.0013;
        params.beta.g = MutualInhibition_beta_g; % old version rho=alpha/alpha: 0.0017;
        
    case 'Neur lateral inhibition zone'
%         params.beta.g = 6*100*0.0013; % 12*100*0.0013;
        params.beta.g = LateralInhibition_beta_g_factor*MutualInhibition_beta_g;
        
    case 'Gaussian-like proneural genes'
        Gaussian_steps = [1 0.61 0.14]; % normalized Gaussian with sigma=1 (1 unit is approximately 1 cell diameter) 
% % %         Gaussian_steps = [1 0.61 0.165];
%         Gaussian_steps = [0.5 0.25 0.025];
        k = params.P*params.Q; % number of cells
        C = getconnectivityM(params.P,params.Q); % which cell connects to which
        SOP_cell = 0.5*(k + params.P); % Choose a cell at around the lattice center
        [~,SOP_nearest_neighbors] = find(C(SOP_cell,:) == 1);
        [~,SOP_next_nearest_neighbors] = find(C(SOP_nearest_neighbors',:) == 1);
        SOP_next_nearest_neighbors = setdiff(SOP_next_nearest_neighbors,[SOP_cell ; SOP_nearest_neighbors']);
%         params.beta.g = 0.0013.*ones(k,1);
%         params.beta.g(SOP_cell) = Gaussian_steps(1).*6*100*0.0013; % 12*100*0.0013; 
%         params.beta.g(SOP_nearest_neighbors) = Gaussian_steps(2).*6*100*0.0013; % 12*100*0.0013;
%         params.beta.g(SOP_next_nearest_neighbors) = Gaussian_steps(3).*6*100*0.0013;
        params.beta.g = MutualInhibition_beta_g.*ones(k,1);
        params.beta.g(SOP_cell) = Gaussian_steps(1).*LateralInhibition_beta_g_factor*MutualInhibition_beta_g; 
        params.beta.g(SOP_nearest_neighbors) = Gaussian_steps(2).*LateralInhibition_beta_g_factor*MutualInhibition_beta_g; 
        params.beta.g(SOP_next_nearest_neighbors) = Gaussian_steps(3).*LateralInhibition_beta_g_factor*MutualInhibition_beta_g;
%         a = 1;
% % %         %%% asymmetric ac\sc expression
% % %         asymetric_exppression = 0.165; %1.18*Gaussian_steps(3);
% % %         a = sort(SOP_next_nearest_neighbors);
% % %         outer_3_left = a(1:3);
% % %         params.beta.g(outer_3_left) = asymetric_exppression.*LateralInhibition_beta_g_factor*MutualInhibition_beta_g;

        params.SOP.cell = SOP_cell;
        params.SOP.neares_neighbors = SOP_nearest_neighbors;
        params.SOP.next_nearest_neighbors = SOP_next_nearest_neighbors;
        
    otherwise
        warning('Unexpected model type. See SOP_defaultparams description')
end

end


%% additional functions
function C=getconnectivityM(P,Q)

k=P*Q; %number of cells
C=zeros(k,k); %This is the connectivity matrix
w=1; %1/6; % Weight for interactions

% calculating the connectivity matrix
for s=1:k
    kneighbour=findneighbourhex(s,P,Q); %finds the neighbors of cell s
    for r=1:6
        C(s,kneighbour(r))=w;
    end
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
end
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