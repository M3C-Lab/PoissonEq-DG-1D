clc;clear all;close all;

% para_Pen_BC = para_IP = 10 / hh

spatial_hh = [0.125000000000000;
              0.062500000000000;
              0.031250000000000;
              0.015625000000000;
              ]; 

p1_H2 = [0.001061114615026;
         2.765219777092091e-04;
         7.051859922631821e-05;
         1.780167347435882e-05;
         ];

p2_H3 = [7.935164505649092e-05;
         1.828848635024315e-05;
         4.386209239383999e-06;
         1.074161609762115e-06;
         ];

log_hh = zeros(length(spatial_hh), 1);
log_p1_H2 = log_hh;
log_p2_H3 = log_hh;

for ii  = 1 : length(spatial_hh)
    log_hh(ii) = log10(spatial_hh(ii));
    log_p1_H2(ii) = log10(p1_H2(ii));
    log_p2_H3(ii) = log10(p2_H3(ii));
end

figure(1)
plot(log_hh, log_p1_H2, "LineWidth", 2);
hold on;
plot(log_hh, log_p2_H3, "LineWidth", 2);
hold on
legend("p1, L2norm / H2norm", "p2, L2norm / H3norm");
title("NIPG");
xlabel("log(hh)");
ylabel("log(relative error)");

a = (log_p1_H2(end) - log_p1_H2(end - 1)) / (log_hh(end) - log_hh(end - 1))
b = (log_p2_H3(end) - log_p2_H3(end - 1)) / (log_hh(end) - log_hh(end - 1))
