function fig5_feature
% Fig. 5 -- Fluctuations in attended feature


% Parameters
neurons = 50;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
beta = 0.1;         % attentional gain
q = 10 / 180 * pi;  % SD[psi] = 10 deg
kappa = 2;
fmean = 10;
f = exp(kappa .* (cos(phi) - 1));   % tuning function
f = f / mean(f) * fmean;
h = cos(phi);
hp = -sin(phi);

% [A] Covariance matrix for small fluctuations around presented direction 
fig = Figure(5, 'size', [200 35]);
subplot(1, 3, 1)
u = q * beta * hp .* f;
C = diag((1 + beta * h) .* f) + (u' * u);
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


% [B] Correlation matrix
subplot(1, 3, 2)
R = C ./ sqrt(diag(C) * diag(C)');
imagesc([-180 180], [-180 180], R)
colormap(cm)
colorbar
caxis([-1, 1 + 2/p] * max(offdiag(R)))
axis square
xlabel('\phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [C] Dependence of correlations on difference in preferred direction for
%     different tuning widths
neurons = 500;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
h = cos(phi);
hp = -sin(phi);
colors = [0.9 0.7 0; 1 0.4 0; 1 0 0; 0.6 0 0.8; 0 0 0.8];
kappa = 2 .^ (-1 : 3);
subplot(1, 3, 3)
hold on
for i = 1 : numel(kappa)
    f = exp(kappa(i) .* (cos(phi) - 1));   % tuning function
    f = f / mean(f) * fmean;
    u = q * beta * hp .* f;
    C = diag((1 + beta * h) .* f) + (u' * u);
    R = C ./ sqrt(diag(C) * diag(C)');
    dphi = abs(angle(exp(1i * bsxfun(@minus, phi, phi')))) / pi * 180;
    dphi = offdiag(dphi);
    r = offdiag(R);
    [ravg, bins] = makeBinned(dphi, r, -0.1 : (360 / neurons) : 180, @nanmean, 'include');
    plot(bins, ravg, 'color', colors(i, :))
end
set(gca, 'XTick', 0 : 90 : 180, 'XLim', [0 180], 'YLim', [-1 1] * 0.002)
xlabel('Difference of preferred directions')
ylabel('Average correlation')
axis square


fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
