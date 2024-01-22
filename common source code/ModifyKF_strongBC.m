function [K, F] = ModifyKF_strongBC(K, F, ID, sol_base)

nNode = length(ID);
for ii = 1 : nNode
    if ID(ii) <= 0
        K(ii, :) = zeros(nNode, 1);
        K(:, ii) = zeros(1, nNode);
        K(ii, ii) = 1.0;
        F(ii) = sol_base(ii);
    end
end

end

