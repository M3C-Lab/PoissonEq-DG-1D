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
nElem = 8;

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

% Symmetric or non-symmetric consistency terms
para_Adj = 1.0;

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

    % Nitsche method on the boundary, assume nElem > 1
    if ee == 1 || ee == nElem
        if ee == 1
            normal_vector = -1.0;
            xi = -1.0;
            g_value = g_left;
        elseif ee == nElem
            normal_vector = 1.0;
            xi = 1.0;
            g_value = g_right;
        end

        k_ele_Nitsche = zeros(nLocBas, nLocBas);
        f_ele_Nitsche = zeros(nLocBas, 1);
        Element = LineElement(ele_order, x_ele, xi);

        for aa = 1 : nLocBas
            NA = Element.Basis(aa);
            NA_x = Element.dN_dx(aa);
            
            f_ele_Nitsche(aa) = f_ele_Nitsche(aa) ...
                - para_Adj * NA_x * normal_vector * g_value ...
                + para_Pen * NA * g_value;

            for bb = 1 : nLocBas
                NB = Element.Basis(bb);
                NB_x = Element.dN_dx(bb);
                
                k_ele_Nitsche(aa, bb) = k_ele_Nitsche(aa, bb) ...
                    - NA * NB_x * normal_vector ...
                    - para_Adj * NA_x * NB * normal_vector ...
                    + para_Pen * NA * NB;
            end
        end
        k_ele = k_ele + k_ele_Nitsche;
        f_ele = f_ele + f_ele_Nitsche;
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