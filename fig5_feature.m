function fig5_feature
% Fig. 5 -- Fluctuations in attended feature


% Parameters
neurons = 50;
phi = linspace(-pi, pi, neurons + 1);
phi = phi(1 : neurons);
beta = 0.1;                 % attentional gain
sdpsi = 10 / 180 * pi;      % variability of attended direction (SD: 10deg)
h = cos(phi);               % gain profile
kappa = 2;                  % tuning width
fmean = 10;                 % mean firing rate across the population
f = exp(kappa .* cos(phi)); % tuning function
f = f / mean(f) * fmean;
mu = exp(beta * h) .* f;    % expected firing rate
mup = kappa * -sin(phi) .* mu;  % derivative of tuning curve

% [A] Covariance matrix for small fluctuations around presented direction 
fig = Figure(5, 'size', [200 35]);
subplot(1, 3, 1)
C = diag(mu) + sdpsi^2 * beta^2 / kappa^2 * (mup' * mup);
imagesc([-180 180], [-180 180], C)
colormap(bluered)
colorbar
caxis([-1, 1] * max(offdiag(C)) * 1.02)
axis square
xlabel('\phi_i')
ylabel('\phi_j')
set(gca, 'XTick', -180 : 90 : 180, 'YTick', -180 : 90 : 180)


% [B] Correlation matrix
subplot(1, 3, 2)
R = C ./ sqrt(diag(C) * diag(C)');
imagesc([-180 180], [-180 180], R)
colormap(bluered)
colorbar
caxis([-1, 1] * max(offdiag(R)) * 1.02)
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
colors = [0.9 0.7 0; 1 0.4 0; 1 0 0; 0.6 0 0.8; 0 0 0.8];
kappa = 2 .^ (-1 : 3);
subplot(1, 3, 3)
hold on
for i = 1 : numel(kappa)
    f = exp(kappa(i) .* cos(phi));   % tuning function
    f = f / mean(f) * fmean;
    mu = exp(beta * h) .* f;
    mup = kappa(i) * -sin(phi) .* mu;
    C = diag(mu) + sdpsi^2 * beta^2 / kappa(i)^2 * (mup' * mup);
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
