function fig3_feature
% Fig. 3 -- Illustration of feature attention model


% Parameters
neurons = 1000;
phi = linspace(-pi, pi, neurons);
deg = linspace(-180, 180, neurons);
nu = 0.2;
kappa = 1;
b = cos(phi);

% attended direction
psi = 0 : pi / 8 : pi;
k = numel(psi);
colors = [1 0 0; 1 0 0.4; 1 0 0.8; 0.8 0 0.9; 0.6 0 1; 0.3 0 1; 0 0 1; 0 0.4 1; 0 0.6 1];

% [A] Individual neuron's tuning function
fig = Figure(3, 'size', [80 70]);
subplot(3, 3, 1)
hold on
f = exp(kappa * (cos(phi) - 1));
for i = 1 : k
    lambda = (1 + nu * cos(psi(i))) * f;
    plot(deg, lambda, 'Color', colors(i, :))
end
plot(deg, f, 'k', 'LineWidth', 2)
set(gca, 'XLim', deg([1 end]), 'XTick', [-180 0 180], 'YTick', [0 1], 'YLim', [0, 1 + nu])
xlabel('Stimulus \theta')
ylabel('\lambda_0(\theta)')
axis square

% [B] Attentional profile b
subplot(3, 3, 2), cla
plot(deg, b, 'k')
set(gca, 'XTick', [-180 0 180], 'XLim', deg([1 end]), 'Ylim', [-1 1])
xlabel('\phi_i - \psi')
ylabel('b(\phi_i - \psi)')
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
    b = cos(phi - psi(i));
    lambda = (1 + nu * b) .* f;
    plot(deg, lambda, 'Color', colors(i, :))
    [m, ndx] = max(lambda);
    plot(deg(ndx), m, '.', 'Color', colors(i, :))
    plot(psi(i) / pi * 180, 0.05, '^', 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :))
end
plot(deg, f, 'k', 'linewidth', 2)
plot([0 0], ylim, ':k')
set(gca, 'XLim', deg([1 end]), 'XTick', -180 : 90 : 180, 'YLim', [0, 1 + nu], 'YTick', [0, 1 - nu, 1, 1 + nu])
xlabel('Preferred direction')
ylabel('\lambda_i(0)')

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
