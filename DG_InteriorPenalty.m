% This code performs well as expected.

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
[node, nNode, hh, ID, IEN] = mesh_1d_DG(Omega_l, Omega_r, nElem, ele_order);

% Quadrature rule
[qp, wq] = Gauss(nLocBas, -1, 1);
nqp = length(qp);

% Penalty parameters
para_Pen_BC = 10 / hh;  % Penalty for weak BC
para_IP = para_Pen_BC;  % Interior penalty

% Consistency parameter
para_Con = 1.0;

% Adjoint-consistency parameter
para_Adj = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIPG: para_IP > some C, para_Con = para_Adj = 1

% NIPG: para_IP > 0, para_Con = 1, para_Adj = -1

% Baumann-Oden: para_IP = 0, para_Con = 1, para_Adj = -1

% Babuska–Zlamal: para_IP > some C, para_Con = para_Adj = 0

% According to the unified formulation, para_IP should be equal to para_Pen_BC,
% However, for the Baumann-Oden method, if para_Pen_BC = 0, it will fail with p1 elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sol_base = zeros(nNode, 1);

% If strongly enforce Dirichlet BC, uncomment the following four lines
% and comment the weak enforcement

sol_base(1) = g_left;
% ID(1) = -ID(1);
% sol_base(end) = g_right;
% ID(end) = -ID(end);

% Initialize the stiffness metrix and load vetor
K = zeros(nNode, nNode);
F = zeros(nNode, 1);

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
                + para_Pen_BC * NA * g_value;

            for bb = 1 : nLocBas
                NB = Element.Basis(bb);
                NB_x = Element.dN_dx(bb);

                k_ele_Nitsche(aa, bb) = k_ele_Nitsche(aa, bb) ...
                    - para_Con * NA * NB_x * normal_vector ...
                    - para_Adj * NA_x * NB * normal_vector ...
                    + para_Pen_BC * NA * NB;
            end
        end
        k_ele = k_ele + k_ele_Nitsche;
        f_ele = f_ele + f_ele_Nitsche;
    end

    % Global assembly
    [K, F] = GAssem(K, F, k_ele, f_ele, ID, IEN(:, ee), sol_base);
    
    % Local assembly for each interface between elements
    if ee ~= nElem
        x_ele_right = zeros(1, nLocBas);
        for aa = 1 : nLocBas
            x_ele_right(aa) = node(IEN(aa, ee + 1));
        end
        
        % Build element basis at the interface
        Element_left = LineElement(ele_order, x_ele, 1.0);
        Element_right = LineElement(ele_order, x_ele_right, -1.0);

        n_left = 1.0;
        n_right = -1.0;

        k_ele_ll = zeros(nLocBas, nLocBas);
        k_ele_rr = zeros(nLocBas, nLocBas);
        k_ele_lr = zeros(nLocBas, nLocBas);
        k_ele_rl = zeros(nLocBas, nLocBas);

        for aa = 1 : nLocBas
            NA_l = Element_left.Basis(aa);
            NA_x_l = Element_left.dN_dx(aa);
            NA_r = Element_right.Basis(aa);
            NA_x_r = Element_right.dN_dx(aa);

            for bb = 1 : nLocBas
                NB_l = Element_left.Basis(bb);
                NB_x_l = Element_left.dN_dx(bb);
                NB_r = Element_right.Basis(bb);
                NB_x_r = Element_right.dN_dx(bb);

                k_ele_ll(aa, bb) = k_ele_ll(aa, bb) ...
                    - 0.5 * n_left * (para_Con * NA_l * NB_x_l + para_Adj * NA_x_l * NB_l) ...
                    + para_IP * n_left * NA_l * n_left * NB_l;

                k_ele_rr(aa, bb) = k_ele_rr(aa, bb) ...
                    - 0.5 * n_right * (para_Con * NA_r * NB_x_r + para_Adj * NA_x_r * NB_r) ...
                    + para_IP * n_right * NA_r * n_right * NB_r;

                k_ele_lr(aa, bb) = k_ele_lr(aa, bb) ...
                    - 0.5 * (para_Con * n_left * NA_l * NB_x_r + para_Adj * NA_x_l * NB_r * n_right) ...
                    + para_IP * n_left * NA_l * n_right * NB_r;

                k_ele_rl(aa, bb) = k_ele_rl(aa, bb) ...
                    - 0.5 * (para_Con * n_right * NA_r * NB_x_l + para_Adj * NA_x_r * NB_l * n_left) ...
                    + para_IP * n_right * NA_r * n_left * NB_l;
            end
        end
        
        % Global assembly on interface
        K = GAssem_interface(K, k_ele_ll, k_ele_rr, k_ele_lr, k_ele_rl, ID, IEN(:, ee), IEN(:, ee + 1), sol_base);
    end
end

[K, F] = ModifyKF_strongBC(K, F, ID, sol_base);

% Solve
uh = K \ F;

% Postprocess
[abs_error, u_normL2, u_normH1, u_normH2, u_normH3] = Postprocess(node, uh, u_exact, u_x, u_xx, u_xxx, Omega_l, Omega_r, nLocBas, nElem, IEN);

rel_error_H2 = abs_error / u_normH2

rel_error_H3 = abs_error / u_normH3

% EOF