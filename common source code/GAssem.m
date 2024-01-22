function [K, F] = GAssem(K, F, k_ele, f_ele, ID, local_IEN, sol_base)
% Global assembly

nLocBas = length(f_ele);
for aa = 1 : nLocBas
    LM_a = ID(local_IEN(aa));
    if LM_a > 0
        F(LM_a) = F(LM_a) + f_ele(aa);

        for bb = 1 : nLocBas
            LM_b = ID(local_IEN(bb));
            if LM_b > 0
                K(LM_a, LM_b) = K(LM_a, LM_b) + k_ele(aa, bb);
            else
                F(LM_a) = F(LM_a) - k_ele(aa, bb) * sol_base(local_IEN(bb));
            end
        end
    end
end

end

