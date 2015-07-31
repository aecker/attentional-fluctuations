function fig7_cm2011
% Fig. 7 -- Modeling Cohen & Maunsell 2011, feature attention

rng(6347563)

% Tuning curve
neurons = 16;
a = 2;
phi = (0 : neurons) / neurons * 2 * pi - pi;
[Phi1, Phi2] = meshgrid(phi);
f = exp(cos(Phi1) + cos(Phi2) + a);

fig = Figure(7, 'size', [150 130]);
subplot(4, 3, 1 : 2)
hdl = surf(Phi1 / pi * 180, Phi2 / pi * 180, f);
set(hdl, 'FaceColor', 'w')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], ...
    'ylim', [-1 1] * 180, 'ytick', [-180 0 180], ...
    'zlim', [0 60], 'ztick', 0 : 20 : 60, 'PlotBoxAspectRatio', [1 1 0.6])
xlabel('Unattended feature')
ylabel('Attended feature')
zlabel('Firing rate (spikes/s)')



% Model I: increased gain variance in unattended condition, but same
%          distribution of attended features
neurons = 50;
a = 2;
phi = (0.5 : neurons) / neurons * 2 * pi - pi;
[Phi1, Phi2] = meshgrid(phi);
f = exp(cos(Phi1(:)) + cos(Phi2(:)) + a);

mu = [0.1 0];
sigma = [0.05 0.1];
stdpsi = 10 / 180 * pi;
k = 25;
psi = linspace(-5, 5, k) * stdpsi;
ppsi = normpdf(psi, 0, stdpsi);
ppsi = ppsi / sum(ppsi);

[~, Cy] = statistics(f, Phi1, Phi2, mu, sigma, psi, ppsi, psi, ppsi);
[Cy1, Cy2] = marginalize(Cy);

subplot(8, 6, 13)
plot(psi / pi * 180, ppsi, 'k')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], 'ytick', [])
xlabel('Attended feature')

subplot(8, 6, 19)
beta = linspace(0.75, 1.25, 50);
plot(beta, normpdf(beta, 1 + mu(1), sigma(1)), 'k')
axis tight
set(gca, 'xlim', beta([1 end]), 'xtick', [0.75 1 1.25], 'ytick', [])
xlabel('Gain')

subplot(8, 6, 25)
plot(psi / pi * 180, ppsi, 'k')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], 'ytick', [])
xlabel('Attended feature')

subplot(8, 6, 31)
beta = linspace(0.75, 1.25, 50);
plot(beta, normpdf(beta, 1 + mu(2), sigma(2)), 'k')
axis tight
set(gca, 'xlim', beta([1 end]), 'xtick', [0.75 1 1.25], 'ytick', [])
xlabel('Gain')

p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
o = ones(p, 1);
cm = [g g o; o gi gi; 0 0 0];
ca = max(abs(offdiag(Cy2))) * (p + 2) / p;

subplot(4, 3, 5)
imagesc([-180 180], [-180 180], Cy1)
caxis([-1 1] * ca)
axis square
colormap(gca, cm)
colorbar
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
xlabel('Neuron i')
ylabel('Neuron j')

subplot(4, 3, 8)
imagesc([-180 180], [-180 180], Cy2)
caxis([-1 1] * ca)
axis square
colormap(gca, cm)
colorbar
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
xlabel('Neuron i')
ylabel('Neuron j')



% Simulated data: using heterogeneous tuning curves
[Rd, binc] = heterogeneous(mu, sigma, psi, ppsi, psi, ppsi);

p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
z = zeros(p, 1);
cmc = [z, 0.75 * gi, 0.75 * gi; g g z];

subplot(4, 3, 9)
imagesc(binc, binc, Rd)
axis square xy
colormap(gca, cmc)
ca2 = 0.12;
caxis([-1 1] * ca2)
colorbar


% Model II: same gain variance in unattended condition, but random feature
%           attended
sigma(2) = sigma(1);
mu(2) = mu(1);
psi2 = (0 : k - 1) / k * 2 * pi;
ppsi2 = ones(1, k) / k;

[~, Cy] = statistics(f, Phi1, Phi2, mu, sigma, psi, ppsi, psi2, ppsi2);
[~, Cy2] = marginalize(Cy);

subplot(8, 6, 37)
plot([-1 1] * 180, [1 1], 'k')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], 'ylim', [0 1.5], 'ytick', [])
xlabel('Attended feature')

subplot(8, 6, 43)
beta = linspace(0.75, 1.25, 50);
plot(beta, normpdf(beta, 1 + mu(2), sigma(2)), 'k')
axis tight
set(gca, 'xlim', beta([1 end]), 'xtick', [0.75 1 1.25], 'ytick', [])
xlabel('Gain')

subplot(4, 3, 11)
imagesc([-180 180], [-180 180], Cy2)
caxis([-1 1] * ca)
axis square
colormap(gca, cm)
colorbar
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
xlabel('Neuron i')
ylabel('Neuron j')


% Simulated data: using heterogeneous tuning curves
[Rd, binc] = heterogeneous(mu, sigma, psi, ppsi, psi2, ppsi2);

subplot(4, 3, 12)
imagesc(binc, binc, Rd)
axis square xy
colormap(gca, cmc)
caxis([-1 1] * ca2)
colorbar
xlabel('Neuron 1 rate change (sp/s)')
ylabel('Neuron 2 rate change (spikes/s)')

fig.cleanup()
fig.save([mfilename('fullpath'), '_.eps'])



function [Ey, Cy] = statistics(f, Phi1, Phi2, mu, sigma, psi1, ppsi1, psi2, ppsi2)

% mean and variance (across attended features) of gain profile
h = @(Phi, psi) cos(bsxfun(@minus, Phi(:), psi));
h1 = h(Phi1, psi1);
Eh1 = h1 * ppsi1';
Eh1h1 = h1 * bsxfun(@times, ppsi1, h1)';
h2 = h(Phi2, psi2);
Eh2 = h2 * ppsi2';
Eh2h2 = h2 * bsxfun(@times, ppsi2, h2)';

% mean and variance of gain
Eg = 1 + (mu(1) * h1) * ppsi1' + (mu(2) * h2) * ppsi2';
Cg = (mu(1)^2 + sigma(1)^2) * Eh1h1 - mu(1)^2 * (Eh1 * Eh1') ...
    + (mu(2)^2 + sigma(2)^2) * Eh2h2 - mu(2)^2 * (Eh2 * Eh2');

% average spike count
Ey = Eg .* f;

% spike count correlations
Cy = diag(Ey) + Cg .* (f * f');
Cy = corrcov(Cy, true);



function [Cy1, Cy2] = marginalize(Cy)

n = sqrt(size(Cy, 1));
Cy1 = 0;
Cy2 = 0;
for j = 1 : n
    for k = 1 : n
        Cy1 = Cy1 + Cy(j : n : end, k : n : end);
        Cy2 = Cy2 + Cy((j - 1) * n + (1 : n), (k - 1) * n + (1 : n));
    end
end
Cy1 = Cy1 / n^2;
Cy2 = Cy2 / n^2;
Cy1(1 : n + 1 : end) = 1;
Cy2(1 : n + 1 : end) = 1;




function [Rd, binc] = heterogeneous(mu, sigma, psi1, ppsi1, psi2, ppsi2)

neurons = 64;
iters = 8;
phi = (0.5 : neurons) / neurons * 2 * pi - pi;
[Phi1, Phi2] = meshgrid(phi);

dmi = zeros(neurons^2 * (neurons^2 - 1) / 2, iters);
dmj = zeros(neurons^2 * (neurons^2 - 1) / 2, iters);
dR = zeros(neurons^2 * (neurons^2 - 1) / 2, iters);
for iter = 1 : iters
    a1 = 0.3 * exp(randn(neurons^2, 1) * 0.3);
    a2 = 0.3 * exp(randn(neurons^2, 1) * 0.3);
    a3 = 3 + randn(neurons^2, 1) * 0.5;
    f = exp(a1 .* cos(Phi1(:)) + a2 .* cos(Phi2(:)) + a3);
    [m1, C] = statistics(f, Phi1, Phi2, mu, sigma, psi1, ppsi1, psi2, ppsi2);
    R1 = corrcov(C, true);
    [m2, C] = statistics(f, Phi1, Phi2, mu([2 1]), sigma([2 1]), psi2, ppsi2, psi1, ppsi1);
    R2 = corrcov(C, true);
    
    [di, dj] = meshgrid(m2 - m1);
    dmi(:, iter) = offdiag(di);
    dmj(:, iter) = offdiag(dj);
    dR(:, iter) = offdiag(R2 - R1);
end

bins = (-5 : 5) * 2;
[Rd, binc] = makeBinned2(dmi(:), dmj(:), dR(:), bins, bins, @mean, 'include');

