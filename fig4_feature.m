function fig4_feature
% Fig. 4 -- Fluctuations in feature attention gain


% Parameters
neurons = 50;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
deg = phi / pi * 180;
nu = 0.1;       % average attentional gain
tau = 0.1;      % variability in attentional gain
kappa = 2;      % tuning width
fmean = 10;
f = exp(kappa .* (cos(phi) - 1));   % tuning function
f = f / mean(f) * fmean;
b = cos(phi);

fig = Figure(4, 'size', [350 35]);
fmax = 100;
xl = [0 fmax];


% [A] Dependence of variance on mean firing rate
subplot(1, 5, 1)
hold on
dp = pi / 4;
psi = 0 : dp : pi;
k = numel(psi);
colors = [0 0 0; 0 0.6 0; 1 0 0; 1 0.7 0; 0 0 1];
for i = 1 : k
    bi = cos(psi(i));
    Ea = linspace(0, fmax, 100);
    Va = Ea + (tau^2 * bi.^2 ./ (1 + bi * nu).^2) * Ea.^2;
    plot(Ea, Va, 'Color', colors(i, :))
end
axis equal
set(gca, 'XLim', xl, 'YLim', [0 200], 'XTick', 0 : 50 : fmax)
xlabel('Mean E[yi]')
ylabel('Variance Var[yi]')

% Illustration of neuron's preferred orientations corresponding to colors
% in panel A
subplot(1, 5, 2)
hold on
plot(deg, f, 'k')
for i = 1 : k
    [~, ndx] = min(abs(psi(i) - phi));
    plot(deg(ndx), f(ndx), '.', 'Color', colors(i, :), 'MarkerSize', 8)
end
set(gca, 'XLim', [-180 180], 'XTick', [-180 0 180], 'YTick', [])
xlabel('Pref. dir. \phi')
ylabel('f_i(0)')


% [B] Covariance structure
subplot(1, 5, 3)
u = tau * b .* f;
C = diag((1 + nu * b) .* f) + (u' * u);
imagesc([-180 180], [-180 180], C)
p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
o = ones(p, 1);
cm = [g g o; o gi gi; 0 0 0];
colormap(cm)
colorbar
caxis([-1, 1 + 2/p] * max(offdiag(C)))
axis square
xlabel('\phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [C] Correlation structure
subplot(1, 5, 4)
R = C ./ sqrt(diag(C) * diag(C)');
imagesc([-180 180], [-180 180], R)
colormap(cm)
colorbar
caxis([-1, 1 + 2/p] * max(offdiag(R)))
axis square
xlabel('\phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [D] Dependence of correlations on difference in preferred direction for
%     different tuning widths
neurons = 500;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
b = cos(phi);
colors = [0.9 0.7 0; 1 0.4 0; 1 0 0; 0.6 0 0.8; 0 0 0.8];
kappa = 2 .^ (-1 : 3);
for i = 1 : numel(kappa)
    f = exp(kappa(i) .* (cos(phi) - 1));   % tuning function
    f = f / mean(f) * fmean;
    u = tau * b .* f;
    C = diag((1 + nu * b) .* f) + (u' * u);
    R = C ./ sqrt(diag(C) * diag(C)');
    dphi = abs(angle(exp(1i * bsxfun(@minus, phi, phi')))) / pi * 180;
    dphi = offdiag(dphi);
    r = offdiag(R);
    [ravg, bins] = makeBinned(dphi, r, -0.1 : (360 / neurons) : 180, @nanmean, 'include');
    
    subplot(1, 5, 5)
    hold on
    plot(bins, ravg, 'color', colors(i, :))
end
subplot(1, 5, 5)
set(gca, 'XTick', 0 : 90 : 180, 'XLim', [0 180], 'YLim', [-1 1] * 0.06)
xlabel('Difference of preferred directions')
ylabel('Average correlation')
axis square


fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
