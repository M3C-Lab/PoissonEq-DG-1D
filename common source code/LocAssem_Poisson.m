function [k_ele, f_ele] = LocAssem_Poisson(nLocBas, Element, f)
% Common k_ele and f_ele for the Poisson equation

k_ele = zeros(nLocBas, nLocBas); f_ele = zeros(nLocBas, 1);

f_qua = f(Element.point_x);

for aa = 1 : nLocBas
    NA = Element.Basis(aa);
    NA_x = Element.dN_dx(aa);

    f_ele(aa) = f_ele(aa) + NA * f_qua;

    for bb = 1 : nLocBas
        NB_x = Element.dN_dx(bb);

        k_ele(aa, bb) = k_ele(aa, bb) + NA_x * NB_x;
    end
end

return;
end

