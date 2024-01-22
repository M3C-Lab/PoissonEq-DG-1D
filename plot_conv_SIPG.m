clc;clear all;close all;

% para_Pen_BC = para_IP = 10 / hh

spatial_hh = [0.125000000000000;
              0.062500000000000;
              0.031250000000000;
              0.015625000000000;
              ]; 

p1_H2 = [0.001318366680871;
         3.353535905729296e-04;
         8.435456783814610e-05;
         2.114088865705560e-05;
         ];

p2_H3 = [6.998177477134878e-06;
         8.696479858826146e-07;
         1.086084330749416e-07;
         1.357766083679906e-08;
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
title("SIPG");
xlabel("log(hh)");
ylabel("log(relative error)");

a = (log_p1_H2(end) - log_p1_H2(end - 1)) / (log_hh(end) - log_hh(end - 1))
b = (log_p2_H3(end) - log_p2_H3(end - 1)) / (log_hh(end) - log_hh(end - 1))
