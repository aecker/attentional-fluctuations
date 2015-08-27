function fig3_feature
% Fig. 3 -- Illustration of feature attention model


% Parameters
neurons = 1000;
phi = linspace(-pi, pi, neurons);
deg = linspace(-180, 180, neurons);
kappa = 1;                  % tuning width
mbeta = 0.2;                % average gain
h = cos(phi);               % gain profile
psi = 0 : pi / 8 : pi;      % attended directions
k = numel(psi);
colors = [1 0 0; 1 0 0.4; 1 0 0.8; 0.8 0 0.9; 0.6 0 1; 0.3 0 1; 0 0 1; 0 0.4 1; 0 0.6 1];


% [A] Individual neuron's tuning function
fig = Figure(3, 'size', [80 70]);
subplot(3, 3, 1)
hold on
f = exp(kappa * (cos(phi) - 1));
for i = 1 : k
    mu = exp(mbeta * cos(psi(i))) * f;
    plot(deg, mu, 'Color', colors(i, :))
end
plot(deg, f, 'k', 'LineWidth', 2)
set(gca, 'XLim', deg([1 end]), 'XTick', [-180 0 180], 'YTick', [0 1], 'YLim', [0, exp(mbeta)])
xlabel('Stimulus \theta')
ylabel('\mu_0(\theta)')
axis square


% [B] Attentional profile b
subplot(3, 3, 2), cla
plot(deg, h, 'k')
set(gca, 'XTick', [-180 0 180], 'XLim', deg([1 end]), 'Ylim', [-1 1])
xlabel('\phi_i - \psi')
ylabel('h(\phi_i - \psi)')
axis square


% [C] Distribution of attentional state term
subplot(3, 3, 3)
x = linspace(-3, 3, 100);
plot(x, normpdf(x), 'k')
set(gca, 'XLim', [-3 3], 'XTick', [-2 0 2])
xlabel('Gain \beta')
ylabel('Density')
axis square


% [D] Population acitivity profile
subplot(3, 3, 4 : 9)
hold all
for i = 1 : k
    h = cos(phi - psi(i));
    mu = exp(mbeta * h) .* f;
    plot(deg, mu, 'Color', colors(i, :))
    [m, ndx] = max(mu);
    plot(deg(ndx), m, '.', 'Color', colors(i, :))
    plot(psi(i) / pi * 180, 0.05, '^', 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :))
end
plot(deg, f, 'k', 'linewidth', 2)
plot([0 0], ylim, ':k')
set(gca, 'XLim', deg([1 end]), 'XTick', -180 : 90 : 180, 'YLim', [0, exp(mbeta)], ...
    'YTick', [0, exp(-mbeta), 1, exp(mbeta)], 'YTickLabel', [0, round(exp(-mbeta) * 10) / 10, 1, round(exp(mbeta) * 10) / 10])
xlabel('Preferred direction \phi_i')
ylabel('Firing rate \mu_i(0)')


fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
