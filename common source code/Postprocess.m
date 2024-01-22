function  [abs_error, u_normL2, u_normH1]  = Postprocess(node, uh, u_exact, u_x, Omega_l, Omega_r, nLocBas, nElem, IEN)
% Postpreocess

% Plot
nPlot_ex = 10001;
x_plot = linspace(Omega_l, Omega_r, nPlot_ex);
u_ex_plot = x_plot;
for ii = 1 : nPlot_ex
    u_ex_plot(ii) = u_exact(x_plot(ii));
end

figure(1);
plot(x_plot, u_ex_plot, "LineWidth", 2);
hold on;
plot(node, uh, "-o", "LineWidth", 2);

xlabel("x-axis");
ylabel("u(x)");
legend("Exact u", "u^h");

% Calculate error
[qp, wq] = Gauss(8, -1, 1);
nqp = length(qp);

abs_error = 0.0;
u_normL2 = 0.0;
u_normH1 = 0.0;

for ee = 1 : nElem
    x_ele = zeros(1, nLocBas);
    uh_ele = zeros(1, nLocBas);
    for aa = 1 : nLocBas
        x_ele(aa) = node((IEN(aa, ee)));
        uh_ele(aa) = uh(IEN(aa, ee));
    end

    for qua = 1 : nqp
        Element = LineElement(nLocBas - 1, x_ele, qp(qua));

        uh_qua = dot(uh_ele, Element.Basis);

        u_qua = u_exact(Element.point_x);
        u_x_qua = u_x(Element.point_x);

        abs_error = abs_error + wq(qua) * Element.Jacobian * (uh_qua - u_qua)^2;
        u_normL2 = u_normL2 + wq(qua) * Element.Jacobian * u_qua^2;
        u_normH1 = u_normH1 + wq(qua) * Element.Jacobian * (u_qua^2 + u_x_qua^2);
    end
end

abs_error = sqrt(abs_error);
u_normL2 = sqrt(u_normL2);
u_normH1 = sqrt(u_normH1);

return;
end

