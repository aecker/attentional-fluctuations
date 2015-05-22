function fig2_spatial
% Fig. 2 -- Spatial attention model


% Parameters
N = 100;
mu = 1.1;                   % average attentional modulation
sigma = [0.05 0.1 0.15];    % variability in attentional gain

fig = Figure(2, 'size', [220 70]);
fmax = 100;
xl = [0 fmax];
colors = [0.4 0 0; 1 0 0; 1 0.6 0.6];

% [A] Dependence of variance on mean firing rate
subplot(2, 4, 1)
f = linspace(0, fmax, N);
hold on
for i = 1 : numel(sigma)
    plot(mu * f, mu * f + sigma(i)^2 * f.^2, 'Color', colors(i, :))
end
plot(xl, xl, ':k')
axis equal
set(gca, 'XLim', xl, 'YLim', [0 200], 'XTick', 0 : 50 : fmax)
xlabel('Mean E[yi]')
ylabel('Variance Var[yi]')

% [B] Dependence of covariances on product of firing rates
subplot(2, 4, 2)
hold on
for i = 1 : numel(sigma)
    plot(mu^2 * f.^2, sigma(i)^2 * f.^2, 'Color', colors(i, :))
end
set(gca, 'XLim', [0 fmax^2], 'Ylim', [0, fmax])
xlabel('Product of means E[yi]E[yj]')
ylabel('Covariance Cov[yi,yj]')
axis square


% [C] Dependence of correlation coefficient on geometric mean firing rate
subplot(2, 4, 3)
hold on
ratio = 2 .^ (0 : 3)';
for i = 1 : numel(sigma)
    for j = 1 : numel(ratio)
        fj = ratio(j) * f;
        rho = sqrt((f .* fj) ./ ((mu / sigma(i)^2 + f) .* (mu / sigma(i)^2 + fj)));
    	hdl(j) = plot(mu * sqrt(f .* fj), rho, 'color', colors(i, :) * (1 - log2(ratio(j)) / numel(ratio))); %#ok
    end
end
set(gca, 'XLim', xl, 'YLim', [0 0.75], 'YTick', 0 : 0.25 : 0.75)
xlabel('Geometric mean rate sqrt(E[yi]E[yj])')
ylabel('Correlation coefficient \rho ij')
axis square
legend(hdl, arrayfun(@(x) sprintf('%d', x), ratio, 'UniformOutput', false))


% [D] Correlation matrix
subplot(2, 4, 5)
N = 50;
phi = linspace(-pi, pi, N + 1); % preferred directions
kappa = 2;                      % tuning width
fmean = 10;                     % average firing rate over all directions
f = exp(kappa * cos(phi));      % tuning function
f = f / mean(f) * fmean;
C = diag(mu * f) + sigma(2)^2 * (f' * f);
R = C ./ sqrt(diag(C) * diag(C)');
imagesc([-180 180], [-180 180], R)
p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
o = ones(p, 1);
cm = [g g o; o gi gi; 0 0 0];
colormap(cm)
colorbar
caxis([-1, 1 + 2/p] * max(offdiag(R)))
axis square
xlabel('\phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [E] Dependence of correlations on difference in preferred direction
subplot(2, 4, 6)
hold on
N = 200;
phi = linspace(-pi, pi, N + 1);
dphi = abs(angle(exp(1i * bsxfun(@minus, phi, phi')))) / pi * 180;
dphi = offdiag(dphi(1 : N, 1 : N));
f = exp(kappa * cos(phi));
f = f / mean(f) * fmean;
bins = (0 : 360 / N : 180) + (180 / N);
for i = 1 : numel(sigma)
    C = diag(mu * f) + sigma(i)^2 * (f' * f);
    R = C ./ sqrt(diag(C) * diag(C)');
    r = offdiag(R(1 : N, 1 : N));
    [ravg, binc] = makeBinned(dphi, r, bins, @nanmean, 'ignore');
    plot(binc, ravg, 'Color', colors(i, :))
end
xlabel('|\phi_i - \phi_j|')
ylabel('Average correlation < \rho ij >')
set(gca, 'XTick', 0 : 45 : 180, 'XLim', [0 180])
axis square


% [F] Dependence of correlations on difference in preferred direction for
%     more narrowly tuned populations
colors = [0.9 0.7 0; 1 0.4 0; 1 0 0; 0.6 0 0.8; 0 0 0.8];
kappa = 2 .^ (-1 : 3);
for i = 1 : numel(kappa)
    f = exp(kappa(i) * cos(phi));
    f = f / mean(f) * fmean;
    C = diag(mu * f) + sigma(2)^2 * (f' * f);
    R = C ./ sqrt(diag(C) * diag(C)');
    dphi = abs(angle(exp(1i * bsxfun(@minus, phi, phi')))) / pi * 180;
    dphi = offdiag(dphi(1 : N, 1 : N));
    r = offdiag(R(1 : N, 1 : N));

    subplot(2, 4, 7)
    hold on
    [ravg, binc] = makeBinned(dphi, r, bins, @nanmean, 'ignore');
    plot(binc, ravg, 'color', colors(i, :))    

    subplot(2, 4, 8)
    hold on
    plot(phi / pi * 180, f, 'color', colors(i, :))
end
subplot(2, 4, 7)
xlabel('|\phi_i - \phi_j|')
ylabel('Average correlation < \rho ij >')
set(gca, 'XTick', 0 : 45 : 180, 'XLim', [0 180], 'YLim', [0 0.15])
axis square
subplot(2, 4, 8)
set(gca, 'XTick', [-180 0 180], 'XLim', [-180 180], 'YLim', [0 max(f)])
plot(xlim, [1 1] * 180, ':k')
axis square

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
