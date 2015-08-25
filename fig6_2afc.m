function fig6_2afc
% Fig. 6 - Correlations due to fluctuating feature attention in 2AFC task

% Parameters
rng(132548)
neurons = 50;
phi = linspace(-pi, pi, neurons + 1);
psi = [-90; 90] / 180 * pi;         % attended directions
beta = 0.1;                         % attentional gain
h = cos(bsxfun(@minus, phi, psi));  % gain profile
f = 10;                             % mean firing rate (sensory response)
mu = exp(beta * h) * f;             % mean rates for attention up/down
dmu = diff(mu);                     % diff between attention conditions


% [A] Illustration of zero-coherence stimulus
fig = Figure(6, 'size', [93 70]);
subplot(2, 2, 1)
n = 40;
x = 2 * rand(n, 2) - 1;
x(sum(x .* x, 2) > 0.9, :) = [];
n = size(x, 1);
ang = rand(n, 1) * 2 * pi;
xx = bsxfun(@plus, x(:, 1), [zeros(n, 1), 0.2 * cos(ang)]);
yy = bsxfun(@plus, x(:, 2), [zeros(n, 1), 0.2 * sin(ang)]);
hold on
plot(x(:, 1), x(:, 2), 'sk', 'MarkerFaceColor', 'k', 'MarkerSize', 2);
plot(xx', yy', 'k');
plot(0.7 * sin(phi), 0.7 * cos(phi), 'k')
axis([-1 1 -1 1])
axis square off


% [B] Sensory and attention-modulated population responses
subplot(2, 2, 2)
colors = [0 0.4 1; 1 0 0];
plot([-180 180], [f f], ':k')
hold on
plot(phi / pi * 180, mu(1, :), 'color', colors(1, :))
plot(phi / pi * 180, mu(2, :), 'color', colors(2, :))
xlabel('Preferred direction')
ylabel('Firing rate (spikes/s)')
legend({'Sensory response', 'Attend left', 'Attend right'}, 'Location', 'SouthWest')
set(gca, 'xlim', [-180 180], 'xtick', -180 : 90 : 180, 'ylim', [0, 1.5 * f])
axis square


% [C] Covariance matrix
subplot(2, 2, 3)
p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
o = ones(p, 1);
cm = [g g o; o gi gi; 0 0 0];
C = f * eye(neurons + 1) + (dmu' * dmu) / 4;
imagesc([-180 180], [-180 180], C)
caxis([-1 1] * 1.1 * max(dmu)^2 / 4)
colormap(cm)
colorbar
xlabel('Preferred direction')
ylabel('Preferred direction')
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
axis square


% [D] Correlation matrix
subplot(2, 2, 4)
R = C ./ sqrt(diag(C) * diag(C)');
imagesc([-180 180], [-180 180], R)
caxis([-1 1] * beta)
colormap(cm)
colorbar
xlabel('Preferred direction')
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
axis square


fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
