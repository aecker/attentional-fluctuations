function fig4_feature
% Fig. 4 -- Fluctuations in feature attention gain


% Parameters
neurons = 50;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
deg = phi / pi * 180;
mbeta = 0.1;                        % average attentional gain
sdbeta = 0.1;                       % variability in attentional gain
h = cos(phi);                       % gain profile
kappa = 2;                          % tuning width
fmean = 10;                         % mean firing rate
f = exp(kappa .* (cos(phi) - 1));   % tuning function
f = f / mean(f) * fmean;

fig = Figure(4, 'size', [350 35]);
fmax = 100;
xl = [0 fmax];


% [A] Dependence of variance on mean firing rate
subplot(1, 5, 1)
hold on
dp = pi / 4;
psi = 0 : dp : pi;
k = numel(psi);
colors = [1 0 0; 1 0 0.8; 0.6 0 1; 0 0 1; 0 0.6 1];
for i = 1 : k
    hi = cos(psi(i));
    mu = linspace(0, exp(mbeta) * fmax, 100);
    plot(mu, mu + sdbeta^2 * (hi .* mu) .^ 2, 'Color', colors(i, :))
end
axis equal
set(gca, 'XLim', xl, 'YLim', [0 200], 'XTick', 0 : 50 : fmax)
xlabel('Mean rate \mu_i')
ylabel('Variance Var[yi]')

% Illustration of neuron's preferred orientations corresponding to colors
% in panel A
subplot(1, 5, 2)
hold on
plot(deg, f, 'k')
for i = 1 : k
    [~, ndx] = min(abs(psi(i) - phi));
    plot(deg(ndx), f(ndx), '.', 'Color', colors(i, :), 'MarkerSize', 15)
end
set(gca, 'XLim', [-180 180], 'XTick', [-180 0 180], 'YTick', [])
xlabel('Pref. dir. \phi')
ylabel('f_i(0)')


% [B] Covariance matrix
subplot(1, 5, 3)
mu = exp(mbeta * h) .* f;
u = h .* mu;
C = diag(mu) + sdbeta^2 * (u' * u);
imagesc([-180 180], [-180 180], C)
colormap(bluered)
colorbar
caxis([-1, 1.02] * max(offdiag(C)))
axis square
xlabel('Preferred direction \phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [C] Correlation matrix
subplot(1, 5, 4)
R = C ./ sqrt(diag(C) * diag(C)');
imagesc([-180 180], [-180 180], R)
colormap(bluered)
colorbar
caxis([-1, 1.02] * max(offdiag(R)))
axis square
xlabel('Preferred direction \phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [D] Dependence of correlations on difference in preferred direction for
%     different tuning widths
neurons = 500;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
h = cos(phi);
colors = [0.9 0.7 0; 1 0.4 0; 1 0 0; 0.6 0 0.8; 0 0 0.8];
kappa = 2 .^ (-1 : 3);
for i = 1 : numel(kappa)
    f = exp(kappa(i) .* (cos(phi) - 1));   % tuning function
    f = f / mean(f) * fmean;
    mu = exp(mbeta * h) .* f;
    u = h .* f;
    C = diag(mu) + sdbeta^2 * (u' * u);
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
