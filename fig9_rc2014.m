function fig9_rc2014()
% Fig. 9 -- Modeling Ruff & Cohen 2014

rng(76348536)

fig = Figure(9, 'size', [150 65]);

% Receptive fields
neurons = 50;
xs = [0.35 0.65];
xc = linspace(0, 1, neurons)';
w = 0.5;
gamma = 3;

subplot(2, 3, 1)
hold on
f = zeros(neurons, 2);
for i = 1 : 2
    f(:, i) = exp(-(xc - xs(i)).^2 ./ (2 * w^2) + gamma);
    plot(xc, f(:, i), ':k')
end
set(gca, 'xtick', [0 xs 1], 'xlim', [0 1], 'ylim', [0 25])
xlabel('Receptive field location')
ylabel('Firing rate (spikes/s)')

% attentional gain profile
x = linspace(-3, 3, 100);
h = @(s, x) 2 * exp(-bsxfun(@minus, x, s) .^ 2 / (2 * w^2)) - exp(-bsxfun(@minus, x, s) .^ 2 / (8 * w^2));

subplot(2, 3, 2)
plot(x, h(0, x), 'k')
set(gca, 'xlim', x([1 end]))
xlabel('Distance from attended location')
ylabel('Gain profile')

% correlation structure
mbeta = 0.2;
sdbeta = 0.05;
sdpsi = [0.5 0.2];  % (unattended, attended)
f = mean(f, 2); % assuming the response to both stimuli simultaneously is
                % the as the average response to the individual stimuli
k = 25;
mu = zeros(size(f));
for i = 1 : 2
    psi = 0.5 + linspace(-5, 5, k) * sdpsi(i);
    ppsi = normpdf(psi, 0.5, sdpsi(i));
    ppsi = ppsi / sum(ppsi);
    [mu(:, i), R] = statistics(f, xc, mbeta, sdbeta, h, psi, ppsi);
    if i == 1
        ca = max(abs(offdiag(R))) * 1.02;
    end
    subplot(2, 3, 6 - i)
    imagesc([0 1], [0 1], R)
    colormap(bluered)
    caxis([-1 1] * ca)
    colorbar
    xlabel('Neuron i')
    ylabel('Neuron j')
    axis square
end

% mean responses: attended and unattended
subplot(2, 3, 1)
col = [0 0 0.7; 0 0.5 0];
plot(xc, mu(:, 1), 'color', col(1, :))
plot(xc, mu(:, 2), 'color', col(2, :))

% focus of attention and distribution of gain
subplot(4, 3, 3)
hold on
xa = linspace(-1, 2, 100);
plot(xa, normpdf(xa, 0.5, sdpsi(1)), 'color', col(1, :))
plot(xa, normpdf(xa, 0.5, sdpsi(2)), 'color', col(2, :))
xlabel('Attended location')
set(gca, 'xlim', [-1 2])

subplot(4, 3, 6)
g = linspace(0, 0.4, 100);
plot(g, normpdf(g, mbeta, sdbeta), 'k')
xlabel('Gain')
set(gca, 'xlim', g([1 end]))

% modeling the experimental data: heterogeneous population
neurons = 512;
xc = linspace(0, 1, neurons)';
w = 0.5 * exp(randn(neurons, 1) * 0.1);
gamma = 2 + randn(neurons, 1) * 0.7;
f = zeros(neurons, 2);
for i = 1 : 2
    f(:, i) = exp(-(xc - xs(i)) .^ 2 ./ (2 * w .^ 2) + gamma);
end
dp = diff(f, [], 2) ./ sqrt(mean(f, 2));
tts = dp * dp';
f = mean(f, 2);

subplot(2, 3, 6)
hold on
bins = -5.5 : 5.5;
for i = 1 : 2
    psi = 0.5 + linspace(-5, 5, k) * sdpsi(i);
    ppsi = normpdf(psi, 0.5, sdpsi(i));
    ppsi = ppsi / sum(ppsi);
    [~, R] = statistics(f, xc, mbeta, sdbeta, h, psi, ppsi);
    [r, binc] = makeBinned(offdiag(tts), offdiag(R), bins, @mean, 'include');
    plot(binc, r, 'color', col(i, :))
end
xlabel('TTS')
ylabel('Correlation')
legend({'Unattended', 'Attended'})
axis square

fig.cleanup()
fig.save([mfilename('fullpath'), '_.eps'])


function [mu, R] = statistics(f, x, mbeta, stbeta, h, psi, ppsi)

% mean and variance (across attended locations) of gain profile
H = h(x, psi);
Eh = H * ppsi';
Ehh = H * bsxfun(@times, ppsi, H)';

% mean and variance of gain
Eg = (mbeta * H) * ppsi';
Cg = (mbeta^2 + stbeta^2) * Ehh - mbeta^2 * (Eh * Eh');

% average spike count
mu = exp(Eg) .* f;

% spike count correlations
C = diag(mu) + Cg .* (mu * mu');
R = C ./ sqrt(diag(C) * diag(C)');

