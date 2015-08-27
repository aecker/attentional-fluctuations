function fig8_cm2011
% Fig. 8 -- Modeling Cohen & Maunsell 2011, feature attention

rng(6347563)

% Illustration of population tuning
neurons = 16;
gamma = 2;
phi = (0 : neurons) / neurons * 2 * pi - pi;
[Phi1, Phi2] = meshgrid(phi);
f = exp(cos(Phi1) + cos(Phi2) + gamma);

fig = Figure(8, 'size', [150 130]);
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
phi = (0.5 : neurons) / neurons * 2 * pi - pi;
[Phi1, Phi2] = meshgrid(phi);
f = exp(cos(Phi1(:)) + cos(Phi2(:)) + gamma);

malpha = [0.1 0];       
sdalpha = [0.05 0.1];
sdpsi = 10 / 180 * pi;
k = 25;
psi = linspace(-5, 5, k) * sdpsi;
ppsi = normpdf(psi, 0, sdpsi);
ppsi = ppsi / sum(ppsi);

[~, R] = statistics(f, Phi1, Phi2, malpha, sdalpha, psi, ppsi, psi, ppsi);
[R1, R2] = marginalize(R);

subplot(8, 6, 13)
plot(psi / pi * 180, ppsi, 'k')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], 'ytick', [])
xlabel('Attended feature')

subplot(8, 6, 19)
beta = linspace(-0.25, 0.25, 50);
plot(beta, normpdf(beta, malpha(1), sdalpha(1)), 'k')
axis tight
set(gca, 'xlim', beta([1 end]), 'xtick', [-0.25 0 0.25], 'ytick', [])
xlabel('Gain')

subplot(8, 6, 25)
plot(psi / pi * 180, ppsi, 'k')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], 'ytick', [])
xlabel('Attended feature')

subplot(8, 6, 31)
plot(beta, normpdf(beta, malpha(2), sdalpha(2)), 'k')
axis tight
set(gca, 'xlim', beta([1 end]), 'xtick', [-0.25 0 0.25], 'ytick', [])
xlabel('Gain')


subplot(4, 3, 5)
imagesc([-180 180], [-180 180], R1)
ca = max(abs(offdiag(R2))) * 1.02;
caxis([-1 1] * ca)
axis square
colormap(gca, bluered)
colorbar
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
xlabel('Neuron i')
ylabel('Neuron j')

subplot(4, 3, 8)
imagesc([-180 180], [-180 180], R2)
caxis([-1 1] * ca)
axis square
colormap(gca, bluered)
colorbar
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
xlabel('Neuron i')
ylabel('Neuron j')


% Simulated data: using heterogeneous tuning curves
[Rd, binc] = heterogeneous(malpha, sdalpha, psi, ppsi, psi, ppsi);

p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
z = zeros(p, 1);
cmc = [z, 0.75 * gi, 0.75 * gi; g g z];

subplot(4, 3, 9)
imagesc(binc, binc, kron(Rd, ones(10)))
axis square xy
colormap(gca, cmc)
ca2 = 0.05;
caxis([-1 1] * ca2)
colorbar


% Model II: same gain variance in unattended condition, but random feature
%           attended
sdalpha(2) = sdalpha(1);
malpha(2) = malpha(1);
psi2 = (0 : k - 1) / k * 2 * pi;
ppsi2 = ones(1, k) / k;

[~, R] = statistics(f, Phi1, Phi2, malpha, sdalpha, psi, ppsi, psi2, ppsi2);
[~, R2] = marginalize(R);

subplot(8, 6, 37)
plot([-1 1] * 180, [1 1], 'k')
set(gca, 'xlim', [-1 1] * 180, 'xtick', [-180 0 180], 'ylim', [0 1.5], 'ytick', [])
xlabel('Attended feature')

subplot(8, 6, 43)
plot(beta, normpdf(beta, malpha(2), sdalpha(2)), 'k')
axis tight
set(gca, 'xlim', beta([1 end]), 'xtick', [-0.25 0 0.25], 'ytick', [])
xlabel('Gain')

subplot(4, 3, 11)
imagesc([-180 180], [-180 180], R2)
caxis([-1 1] * ca)
axis square
colormap(gca, bluered)
colorbar
set(gca, 'xtick', -180 : 90 : 180, 'ytick', -180 : 90 : 180)
xlabel('Neuron i')
ylabel('Neuron j')


% Simulated data: using heterogeneous tuning curves
[Rd, binc] = heterogeneous(malpha, sdalpha, psi, ppsi, psi2, ppsi2);

subplot(4, 3, 12)
imagesc(binc, binc, kron(Rd, ones(10)))
axis square xy
colormap(gca, cmc)
caxis([-1 1] * ca2)
colorbar
xlabel('Neuron 1 rate change (sp/s)')
ylabel('Neuron 2 rate change (spikes/s)')

fig.cleanup()
fig.save([mfilename('fullpath'), '_.eps'])



function [mu, R] = statistics(f, Phi1, Phi2, malpha, sdalpha, psi1, ppsi1, psi2, ppsi2)

% mean and variance (across attended features) of gain profile
h = @(Phi, psi) cos(bsxfun(@minus, Phi(:), psi));
h1 = h(Phi1, psi1);
Eh1 = h1 * ppsi1';
Eh1h1 = h1 * bsxfun(@times, ppsi1, h1)';
h2 = h(Phi2, psi2);
Eh2 = h2 * ppsi2';
Eh2h2 = h2 * bsxfun(@times, ppsi2, h2)';

% mean and variance of attentional gain [g = alpha * h(psi)]
Eg = (malpha(1) * h1) * ppsi1' + (malpha(2) * h2) * ppsi2';
Cg = (malpha(1)^2 + sdalpha(1)^2) * Eh1h1 - malpha(1)^2 * (Eh1 * Eh1') ...
   + (malpha(2)^2 + sdalpha(2)^2) * Eh2h2 - malpha(2)^2 * (Eh2 * Eh2');

% average spike count
mu = exp(Eg) .* f;

% spike count correlations
C = diag(mu) + Cg .* (mu * mu');
R = C ./ sqrt(diag(C) * diag(C)');



function [R1, R2] = marginalize(R)

n = sqrt(size(R, 1));
R1 = 0;
R2 = 0;
for j = 1 : n
    for k = 1 : n
        R1 = R1 + R(j : n : end, k : n : end);
        R2 = R2 + R((j - 1) * n + (1 : n), (k - 1) * n + (1 : n));
    end
end
R1 = R1 / n^2;
R2 = R2 / n^2;
R1(1 : n + 1 : end) = 1;
R2(1 : n + 1 : end) = 1;




function [Rd, binc] = heterogeneous(malpha, sdalpha, psi1, ppsi1, psi2, ppsi2)

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
    [m1, R1] = statistics(f, Phi1, Phi2, malpha, sdalpha, psi1, ppsi1, psi2, ppsi2);
    [m2, R2] = statistics(f, Phi1, Phi2, malpha([2 1]), sdalpha([2 1]), psi2, ppsi2, psi1, ppsi1);
    
    [di, dj] = meshgrid(m2 - m1);
    dmi(:, iter) = offdiag(di);
    dmj(:, iter) = offdiag(dj);
    dR(:, iter) = offdiag(R2 - R1);
end

bins = (-5 : 5) * 2;
[Rd, binc] = makeBinned2(dmi(:), dmj(:), dR(:), bins, bins, @mean, 'include');

