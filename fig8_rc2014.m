function fig8_rc2014()
% Fig. 8 -- Modeling Ruff & Cohen 2014

rng(76348536)

fig = Figure(8, 'size', [150 65]);

% Receptive fields
neurons = 50;
xs = [0.35 0.65];
xc = linspace(0, 1, neurons)';
w = 0.5;
a = 3;

subplot(2, 3, 1)
hold on
f = zeros(neurons, 2);
for i = 1 : 2
    f(:, i) = exp(-(xc - xs(i)).^2 ./ (2 * w^2) + a);
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
mu = 0.2;
sigma = 0.05;
stdpsi = [0.5 0.2];  % (unattended, attended)
f = mean(f, 2); % assuming the response to both stimuli simultaneously is
                % the as the average response to the individual stimuli
p = 100;
g = (1/p : 1/p : 1)';
gi = flipud(g);
o = ones(p, 1);
cm = [g g o; o gi gi; 0 0 0];
k = 25;
Ey = zeros(size(f));
for i = 1 : 2
    psi = 0.5 + linspace(-5, 5, k) * stdpsi(i);
    ppsi = normpdf(psi, 0.5, stdpsi(i));
    ppsi = ppsi / sum(ppsi);
    [Ey(:, i), Cy] = statistics(f, xc, mu, sigma, h, psi, ppsi);
    Cy = corrcov(Cy, true);
    if i == 1
        ca = max(abs(offdiag(Cy))) * (p + 2) / p;
    end
    subplot(2, 3, 6 - i)
    imagesc([0 1], [0 1], Cy)
    colormap(cm)
    caxis([-1 1] * ca)
    colorbar
    xlabel('Neuron i')
    ylabel('Neuron j')
    axis square
end

% mean responses: attended and unattended
subplot(2, 3, 1)
col = [0 0 0.7; 0 0.5 0];
plot(xc, Ey(:, 1), 'color', col(1, :))
plot(xc, Ey(:, 2), 'color', col(2, :))


% modeling the experimental data: heterogeneous population
neurons = 512;
xc = linspace(0, 1, neurons)';
w = 0.5 * exp(randn(neurons, 1) * 0.1);
a = 2 + randn(neurons, 1) * 0.7;
f = zeros(neurons, 2);
for i = 1 : 2
    f(:, i) = exp(-(xc - xs(i)) .^ 2 ./ (2 * w .^ 2) + a);
end
dp = diff(f, [], 2) ./ sqrt(mean(f, 2));
tts = dp * dp';
f = mean(f, 2);

subplot(2, 3, 6)
hold on
bins = -5.5 : 5.5;
R = zeros(neurons, neurons, 2);
for i = 1 : 2
    psi = 0.5 + linspace(-5, 5, k) * stdpsi(i);
    ppsi = normpdf(psi, 0.5, stdpsi(i));
    ppsi = ppsi / sum(ppsi);
    [~, Cy] = statistics(f, xc, mu, sigma, h, psi, ppsi);
    R(:, :, i) = corrcov(Cy, true);
    [r, binc] = makeBinned(offdiag(tts), offdiag(R(:, :, i)), bins, @mean, 'include');
    plot(binc, r, 'color', col(i, :))
end
xlabel('TTS')
ylabel('Correlation')
legend({'Unattended', 'Attended'})
axis square

fig.cleanup()
fig.save([mfilename('fullpath'), '_.eps'])


function [Ey, Cy] = statistics(f, x, mu, sigma, h, psi, ppsi)

% mean and variance (across attended locations) of gain profile
H = h(x, psi);
Eh = H * ppsi';
Ehh = H * bsxfun(@times, ppsi, H)';

% mean and variance of gain
Eg = 1 + (mu * H) * ppsi';
Cg = (mu^2 + sigma^2) * Ehh - mu^2 * (Eh * Eh');

% average spike count
Ey = Eg .* f;

% spike count correlations
Cy = diag(Ey) + Cg .* (f * f');
Cy = corrcov(Cy, true);



