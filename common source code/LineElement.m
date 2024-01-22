function LE = LineElement(ele_order, x_ele, xi)
% Generate the basis function N, dN/dx, Jacobian and physical coordinate in line element with given xi

% Input: x_ele = [x1, x2, x3, ...] of local nodes
%        xi in [-1, 1]

nlocbas = ele_order + 1;

LE.Basis = zeros(nlocbas, 1);

dN_dxi = zeros(nlocbas, 1);
for aa = 1 : nlocbas
    % [N1; N2; N3; ...]
    LE.Basis(aa) = PolyBasis(ele_order, aa, 0, xi);

    % [N1_xi; N2_xi; N3_xi; ...]
    dN_dxi(aa) = PolyBasis(ele_order, aa, 1, xi);
end

% dx_dxi
LE.Jacobian = dot(x_ele, dN_dxi);

% [N1_x; N2_x; N3_x; ...]
LE.dN_dx = dN_dxi / LE.Jacobian;

% x_qua
LE.point_x = dot(x_ele, LE.Basis);

return;
end

