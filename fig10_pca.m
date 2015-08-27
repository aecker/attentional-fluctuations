function fig10_pca
% Fig. 10 -- Subspace spanned by attentional fluctuations


% Parameters
neurons = 50;
phi = linspace(-pi, pi, neurons + 1);
deg = phi / pi * 180;
kappa = 1;
beta = 0.2;
psi = 0;
normalize = @(x) x / norm(x(1 : end - 1));


% [A] Regular principal components
fig = Figure(10, 'size', [100 30]);
f = @(theta) exp(kappa * cos(theta - phi));
fp = @(theta) -kappa * sin(theta - phi) .* f(theta);
h = @(psi) beta * cos(psi - phi);
hp = @(psi) -beta * sin(psi - phi);
i = 1;
for theta = [0, pi / 3]
    subplot(1, 3, i)
    hold on
    cb = [0 0 1; 0 0.4 1];
    plot(deg, normalize(f(theta)) - 0.005, 'k')
    plot(deg, normalize(fp(theta)), 'color', 0.5 * [1 1 1])
    plot(deg, normalize(f(theta)), 'r')
    plot(deg, normalize(h(psi) .* f(theta)), 'color', cb(1, :))
    plot(deg, normalize(hp(psi) .* f(theta)), 'color', cb(2, :))
    set(gca, 'xlim', [-1 1] * 180, 'xtick', (-1 : 0.5 : 1) * 180, ...
        'ylim', [-1 1] * 0.35, 'ytick', 0)
    i = i + 1;
    xlabel('Preferred direction')
    ylabel('Weight')
    axis square
end
legend({'f', 'f''', '\alpha', '\beta', '\psi'})


% [B] Exponential Family PCA
subplot(1, 3, 3)
hold on
plot([-1 1] * 180, [1 1] / sqrt(neurons), 'r')
plot(deg, normalize(cos(psi - phi)), 'color', cb(1, :))
plot(deg, normalize(-sin(psi - phi)), 'color', cb(2, :))
set(gca, 'xlim', [-1 1] * 180, 'xtick', (-1 : 0.5 : 1) * 180, ...
    'ylim', [-1 1] * 0.35, 'ytick', 0)
xlabel('Preferred direction')
axis square


fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
