function [node, nNode, hh, ID, IEN] = mesh_1d_DG(Omega_l, Omega_r, nElem, ele_order)
% Get the uniform node coordinate and mesh size in 1-D problem with DG.
% Generate ID and IEN array, IEN = LM if dof = 1.

if Omega_l >= Omega_r || nElem <= 0 || ele_order <= 0
    disp("mesh_1d_DG: Wrong input.");
    return;
end

hh = (Omega_r - Omega_l) / nElem;

nNode = (ele_order + 1) * nElem;

x_coor = linspace(Omega_l, Omega_r, (ele_order * nElem) + 1);
node = zeros(nNode, 1);

for ee = 1 : nElem
    for aa = 1 : ele_order + 1
        node((ele_order + 1) * (ee - 1) + aa) = x_coor(ele_order * (ee - 1) + aa); 
    end
end

ID = linspace(1, nNode, nNode);
IEN = zeros(ele_order + 1, nElem);

nn = 1;
for ee = 1 : nElem
    for aa = 1 : ele_order + 1
        IEN(aa, ee) = nn;
        nn = nn + 1;
    end
end

return;
end
