function plotFCF(fcf)

figure
plot(fcf(:, 1), fcf(:, 2))
title('Feature Counting Function')
xlabel('epsilon')
ylabel('rank')