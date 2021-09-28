function y0 = SOP_InitialConditions_rho(model,params,k)
%SOP_INITIALCONDITIONS Summary of this function goes here
%   Detailed explanation goes here
        
    noise = 1 + 0.05.*2.*(rand(k,1) - 0.5); % noise up to 5 percent
    d0 = 1.*1.*ones(k,1);
    d_M0 = 1.*0.1.*ones(k,1); % old version rho=alpha/alpha: 0.062.*ones(k,1);
    d_N0 = zeros(k,1);
    n0 = 1.*1.*ones(k,1); % All initial Notch levels are 1
%     E0 = 0.0058.*ones(k,1);
    E0 = 0.01.*ones(k,1); % old version rho=alpha/alpha: 0.0037.*ones(k,1);
%     g0 = model dependent
    MutualInhibition_g = 0.00051;
    LateralInhibition_g_factor = 100;
    N0 = zeros(k,1);
    M0 = ones(k,1);
    switch model
        case 'Mib1 mutual inhibition zone'
%             g0 = 0.00014.*ones(k,1).*noise;
            g0 = MutualInhibition_g.*ones(k,1).*noise;% old version rho=alpha/alpha: 0.00019.*ones(k,1).*noise;

        case 'Neur lateral inhibition zone'
%             g0 = 0.014.*ones(k,1).*noise;
            g0 = LateralInhibition_g_factor.*MutualInhibition_g.*ones(k,1).*noise;

        case 'Gaussian-like proneural genes'
            Gaussian_steps = [1 0.61 0.14];
% % %             Gaussian_steps = [1 0.61 0.165];
            C = getconnectivityM(params.P,params.Q); % which cell connects to which
            SOP_cell = 0.5*(k + params.P); % Choose a cell at around the lattice center
            [~,SOP_nearest_neighbors] = find(C(SOP_cell,:) == 1);
            [~,SOP_next_nearest_neighbors] = find(C(SOP_nearest_neighbors',:) == 1);
            SOP_next_nearest_neighbors = setdiff(SOP_next_nearest_neighbors,[SOP_cell ; SOP_nearest_neighbors']);
%             g0 = 0.00014.*ones(k,1).*noise;
%             g0(SOP_cell) = Gaussian_steps(1).*0.014.*noise(SOP_cell,1);
%             g0(SOP_nearest_neighbors) = Gaussian_steps(2).*0.014.*noise(SOP_nearest_neighbors,1);
%             g0(SOP_next_nearest_neighbors) = Gaussian_steps(3).*0.014.*noise(SOP_next_nearest_neighbors,1);
            g0 = MutualInhibition_g.*ones(k,1).*noise;
            g0(SOP_cell) = Gaussian_steps(1).*LateralInhibition_g_factor.*MutualInhibition_g.*noise(SOP_cell,1);
            g0(SOP_nearest_neighbors) = Gaussian_steps(2).*LateralInhibition_g_factor.*MutualInhibition_g.*noise(SOP_nearest_neighbors,1);
            g0(SOP_next_nearest_neighbors) = Gaussian_steps(3).*LateralInhibition_g_factor.*MutualInhibition_g.*noise(SOP_next_nearest_neighbors,1);
%             a = 1;
% % %             %%% asymmetric ac\sc expression
% % %             asymetric_exppression = 0.165; %1.18*Gaussian_steps(3);
% % %             a = sort(SOP_next_nearest_neighbors);
% % %             outer_3_left = a(1:3);
% % %             g0(outer_3_left) = asymetric_exppression.*LateralInhibition_g_factor.*MutualInhibition_g.*noise(outer_3_left,1);

        otherwise
            warning('Unexpected model type. See SOP_defaultparams description')
    end

    y0 = [d0;d_M0;d_N0;n0;E0;g0;N0;M0]; % vector of initial conditions
    
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