function K = GAssem_interface(K, k_11, k_22, k_12, k_21, ID, local_IEN_1, local_IEN_2, sol_base)
% Global assembly for some point of interface of two elements

nLocBas = length(local_IEN_1);

for aa = 1 : nLocBas
    LM_a1 = ID(local_IEN_1(aa));
    if LM_a1 > 0
        for bb = 1 : nLocBas
            LM_b1 = ID(local_IEN_1(bb));
            if LM_b1 > 0
                K(LM_a1, LM_b1) = K(LM_a1, LM_b1) + k_11(aa, bb);
            else
                K(LM_a1, local_IEN_1(bb)) = K(LM_a1, local_IEN_1(bb)) - k_11(aa, bb) * sol_base(local_IEN_1(bb));
            end

            LM_b2 = ID(local_IEN_2(bb));
            if LM_b2 > 0
                K(LM_a1, LM_b2) = K(LM_a1, LM_b2) + k_12(aa, bb);
            else
                K(LM_a1, local_IEN_2(bb)) = K(LM_a1, local_IEN_2(bb)) - k_12(aa, bb) * sol_base(local_IEN_2(bb));
            end
        end
    end

    LM_a2 = ID(local_IEN_2(aa));
    if LM_a2 > 0
        for bb = 1 : nLocBas
            LM_b1 = ID(local_IEN_1(bb));
            if LM_b1 > 0
                K(LM_a2, LM_b1) = K(LM_a2, LM_b1) + k_21(aa, bb);
            else
                K(LM_a2, local_IEN_1(bb)) = K(LM_a2, local_IEN_1(bb)) - k_21(aa, bb) * sol_base(local_IEN_1(bb));
            end

            LM_b2 = ID(local_IEN_2(bb));
            if LM_b2 > 0
                K(LM_a2, LM_b2) = K(LM_a2, LM_b2) + k_22(aa, bb);
            else
                K(LM_a2, local_IEN_2(bb)) = K(LM_a2, local_IEN_2(bb)) - k_22(aa, bb) * sol_base(local_IEN_2(bb));
            end
        end
    end
end

return;
end

