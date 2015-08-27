function fig2_spatial
% Fig. 2 -- Spatial attention model


% Parameters
N = 100;
malpha = 0.1;               % average attentional modulation
sigma = [0.05 0.1 0.15];    % variability in attentional gain

fig = Figure(2, 'size', [220 70]);
fmax = 100;
xl = [0 fmax];
colors = [0.4 0 0; 1 0 0; 1 0.6 0.6];


% [A] Dependence of variance on mean firing rate
subplot(2, 4, 1)
f = linspace(0, fmax, N);
mu = exp(malpha) * f;
hold on
for i = 1 : numel(sigma)
    plot(mu, mu + sigma(i)^2 * mu .^ 2, 'Color', colors(i, :))
end
plot(xl, xl, ':k')
axis equal
set(gca, 'XLim', xl, 'YLim', [0 200], 'XTick', 0 : 50 : fmax)
xlabel('Mean rate \mu_i')
ylabel('Variance')


% [B] Dependence of covariances on product of firing rates
subplot(2, 4, 2)
hold on
for i = 1 : numel(sigma)
    plot(mu .^ 2, sigma(i)^2 * mu .^ 2, 'Color', colors(i, :))
end
set(gca, 'XLim', [0 fmax^2], 'Ylim', [0, fmax])
xlabel('Product of mean rates \mu_i\mu_j')
ylabel('Covariance Cov[yi,yj]')
axis square


% [D] Dependence of correlation coefficient on geometric mean firing rate
subplot(2, 4, 5)
hold on
ratio = 2 .^ (0 : 3)';
for i = 1 : numel(sigma)
    for j = 1 : numel(ratio)
        muj = ratio(j) * mu;
        rho = sqrt(mu ./ (sigma(i)^-2 + mu) .* muj ./ (sigma(i)^-2 + muj));
    	hdl(j) = plot(sqrt(mu .* muj), rho, 'color', colors(i, :) * (1 - log2(ratio(j)) / numel(ratio))); %#ok
    end
end
set(gca, 'XLim', xl, 'YLim', [0 0.75], 'YTick', 0 : 0.25 : 0.75)
xlabel('Geometric mean rate (\mu_i\mu_j)^{1/2}')
ylabel('Correlation coefficient \rho_{ij}')
axis square
legend(hdl, arrayfun(@(x) sprintf('%d', x), ratio, 'UniformOutput', false))


% [C] Covariance matrix
subplot(2, 4, 3)
N = 50;
phi = linspace(-pi, pi, N + 1); % preferred directions
kappa = 2;                      % tuning width
fmean = 10;                     % average firing rate over all directions
f = exp(kappa * cos(phi));      % tuning function
f = f / mean(f) * fmean;
mu = exp(malpha) .* f;
C = diag(mu) + sigma(2)^2 * (mu' * mu);
imagesc([-180 180], [-180 180], C)
colormap(bluered)
colorbar
caxis([-1, 1.02] * max(offdiag(C)))
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
mu = exp(malpha) * f;
bins = (0 : 360 / N : 180) + (180 / N);
for i = 1 : numel(sigma)
    C = diag(mu) + sigma(i)^2 * (mu' * mu);
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
    mu = exp(malpha) * f;
    C = diag(mu) + sigma(2)^2 * (mu' * mu);
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
