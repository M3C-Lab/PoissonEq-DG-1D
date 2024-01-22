clc; clear all; close all;
addpath("common source code/")

% Domain
Omega_l = 0.0;
Omega_r = 1.0;

% Exact solution
u_exact = @(x) 5 * sin(4 * x) - 3 * x.^3 + 1;
u_x = @(x) 20 * cos(4 * x) - 9 * x.^2;
f = @(x) 80 * sin(4 * x) + 18 * x;
u_xx = @(x) -f(x);
u_xxx = @(x) -320 * cos(4 * x) + 18;

% Dirichlet BC at Omega_l and Omega_r
g_left = u_exact(Omega_l);
g_right = u_exact(Omega_r);

% Number of elements
nElem = 10;

% Element order
ele_order = 1;

% Number of local basis function
nLocBas = ele_order + 1;

% Node coordinate and mesh size, ID and IEN array
[node, nNode, hh, ID, IEN] = mesh_1d_CG(Omega_l, Omega_r, nElem, ele_order);

% Quadrature rule
[qp, wq] = Gauss(nLocBas, -1, 1);
nqp = length(qp);

% Penalty parameter
para_Pen = 10 / hh;

% Initialize the stiffness metrix and load vetor
K = zeros(nNode, nNode);
F = zeros(nNode, 1);

sol_base = zeros(nNode, 1);

% If strongly enforce Dirichlet BC, uncomment the following four lines
% and comment the weak enforcement

% sol_base(1) = g_left;
% ID(1) = -ID(1);
% sol_base(end) = g_right;
% ID(end) = -ID(end);

for ee = 1 : nElem
    % Get local node coordinate
    x_ele = zeros(1, nLocBas);
    for aa = 1 : nLocBas
        x_ele(aa) = node(IEN(aa, ee));
    end

    k_ele = zeros(nLocBas, nLocBas);
    f_ele = zeros(nLocBas, 1);
    
    % Local assembly in the domain
    for qua = 1 : nqp
        Element = LineElement(ele_order, x_ele, qp(qua));
        
        % Common part of the weak form
        [k_ele_Poisson, f_ele_Poisson] = LocAssem_Poisson(nLocBas, Element, f);

        k_ele = k_ele + wq(qua) * Element.Jacobian * k_ele_Poisson;
        f_ele = f_ele + wq(qua) * Element.Jacobian * f_ele_Poisson;
    end

    % Penalty on the boundary, assume nElem > 1
    if ee == 1
        Element = LineElement(ele_order, x_ele, -1.0);
        k_ele(1, 1)  = k_ele(1, 1) + para_Pen * Element.Basis(1) * Element.Basis(1);
        f_ele(1) = f_ele(1) + para_Pen * Element.Basis(1) * g_left;
    elseif ee == nElem
        Element = LineElement(ele_order, x_ele, 1.0);
        k_ele(nLocBas, nLocBas) = k_ele(nLocBas, nLocBas) + para_Pen * Element.Basis(nLocBas) * Element.Basis(nLocBas);
        f_ele(nLocBas) = f_ele(nLocBas) + para_Pen * Element.Basis(nLocBas) * g_right;
    end

    % Global assembly
    [K, F] = GAssem(K, F, k_ele, f_ele, ID, IEN(:, ee), sol_base);
end

[K, F] = ModifyKF_strongBC(K, F, ID, sol_base);

% Solve
uh = K \ F;

% Postprocess
[abs_error, u_normL2, u_normH1, u_normH2, u_normH3] = Postprocess(node, uh, u_exact, u_x, u_xx, u_xxx, Omega_l, Omega_r, nLocBas, nElem, IEN);

rel_error = abs_error / u_normL2

% EOF
