function [node, nNode, hh, ID, IEN] = mesh_1d_CG(Omega_l, Omega_r, nElem, ele_order)
% Get the uniform node coordinate and mesh size in 1-D problem
% Generate ID and IEN array, IEN = LM if dof = 1.

% Note that nNode = nEquation here.
% If strongly enforce Dirichlet BC g_{A} on node A, 
% let F(LM_a) = F(LM_a) - k_ele(aa, bb) * sol_base(LM_b) in GAssem,
% then let K_{AA} = 1, K_{AB} = K{BA} = 0 for B != A and F_{A} = g_{A} in ModifyKF_strongBC.

if Omega_l >= Omega_r || nElem <= 0 || ele_order <= 0
    disp("mesh_1d_CG: Wrong input.");
    return;
end

hh = (Omega_r - Omega_l) / nElem;

nNode = (ele_order * nElem) + 1;

node = linspace(Omega_l, Omega_r, nNode);
node = node';

ID = linspace(1, nNode, nNode);
IEN = zeros(ele_order + 1, nElem);

nn = 1;
for ee = 1 : nElem
    for aa = 1 : ele_order + 1
        IEN(aa, ee) = nn;
        if aa ~= ele_order + 1
            nn = nn + 1;
        end
    end
end

return;
end

