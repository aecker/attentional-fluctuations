function fig1_spatial
% Fig. 1 -- Spatial attention model


% Parameters
N = 100;
phi = linspace(-pi, pi, N);
deg = linspace(-180, 180, N);
mu = 0.2;       % average attentional gain
kappa = 1;      % tuning width

% [A] illustrating an example stimulus
fig = Figure(1, 'size', [80 70]);
subplot(3, 3, 1)
n = 80;
x = 2 * rand(n, 2) - 1;
x(sum(x .* x, 2) > 1, :) = [];
hold on
plot(x(:, 1), x(:, 2), 'sk', 'MarkerFaceColor', 'k', 'MarkerSize', 2);
plot(0.7 * sin(phi), 0.7 * cos(phi), 'k')
axis([-1 1 -1 1])
axis square off

% [B] Individual neuron's tuning function
subplot(3, 3, 2)
f = exp(kappa * (cos(phi) - 1));
plot(deg, f, ':k', deg, (1 + mu) * f, 'k')
set(gca, 'XLim', [-180 180], 'XTick', [-180 0 180], 'YTick', [], 'YLim', [0 (1+mu+0.1)*max(f)], 'Box', 'off')
xlabel('Stimulus \theta')
ylabel('f(\theta)')
axis square

% [C] Distribution of attentional state term
subplot(3, 3, 3)
x = linspace(-3, 3, 100);
plot(x, normpdf(x), 'k')
set(gca, 'XLim', [-3 3], 'XTick', [-2 0 2], 'Box', 'off')
xlabel('Gain \alpha')
ylabel('Density')
axis square

% [D] Population acitivity profile
subplot(3, 3, 4:9)
plot(deg, f, ':k', deg, (1 + mu) * f, 'k')
set(gca, 'XLim', [-180 180], 'XTick', -180 : 90 : 180, 'YTick', [0, 1, 1+mu], 'YLim', [0, 1+mu], 'Box', 'off')
xlabel('Preferred direction')
ylabel('f(\theta)')

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
