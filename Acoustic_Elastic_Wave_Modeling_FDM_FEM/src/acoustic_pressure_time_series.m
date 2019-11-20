plot(t(1: (end-1)), u55(:, 45), 'linewidth', 3, 'color', 'k');
hold on;
plot(t(1: (end-1)), u30(:, 45), 'linewidth', 3, 'color', 'r');
% scatter(t(1: (end-1)), u30(:, 45), 60, '+', ...
%     'linewidth', 1, 'markeredgecolor', 'r');
xlabel('Time [s]');
ylabel('Amp.');
xlim([t(end)*0.4, t(end)]);
legend('55 Hz', '30 Hz');
grid on;
set(gca, 'fontsize', 20, 'fontweight', 'bold');