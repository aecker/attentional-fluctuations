function fig7_fisherinfo
% Fig. 7 - Fisher information under fluctuations of attentional state


% set input noise level for given threshold
criterion = 0.75;               % threshold criterion (% correct)
dp = 2 * norminv(criterion);    % equivalent d'
thresh = [0.5 3] / 180 * pi;    % threshold based on input noise
sdtheta = thresh / dp;          % SD of input noise

n = 10 : 10 : 500;              % number of neurons
malpha = 0.1;                   % average spatial gain
sdalpha = 0.1;                  % variance of spatial gain
mbeta = 0.2;                    % average feature gain
sdbeta = 0.1;                   % variance of feature gain
sdpsi = 10 / 180 * pi;          % variance of attended direction

% tuning parameters
kappa = 2;                      % tuning width
gamma = log(30) - kappa;        % average peak firing rate
amplvar = 1.5;                  % variance of amplitudes across the population

nk = numel(n);
nt = numel(sdtheta);
J = zeros(nk, 4, nt, 2);

% generate population
for l = 1 : 2
    for k = 1 : nk
        neurons = n(k);
        phi = (0.5 : neurons) / neurons * 2 * pi - pi;
        f = exp(kappa .* cos(phi) + gamma);
        df = -kappa .* sin(-phi) .* f;
        h = cos(phi);
        
        if l == 1
            mu = f;
            dmu = df;
        else
            mu = exp(malpha + mbeta * h) .* f;
            dmu = exp(malpha + mbeta * h) .* df;
        end
        
        % (1) no attentional fluctuations
        J0 = dmu * diag(1 ./ mu) * dmu';
        J(k, 1, 1, l) = J0;
        
        % (2) spatial gain fluctuations
        J(k, 2, 1, l) = J0 - amplvar * sum(dmu .^ 2) / (sdalpha ^ -2 + sum(mu));
        
        % (3) feature gain fluctuations
        J(k, 3, 1, l) = J0 - amplvar * sum((h .* dmu) .^ 2) / (sdbeta ^ -2 + sum(h .^ 2 .* mu));
        
        % (4) attended feature fluctuations
        J(k, 4, 1, l) = J0 / (1 + (sdpsi ^ 2 * mbeta ^ 2 / kappa ^ 2) * J0);
        
        % add effect pf input noise
        for j = 1 : nt
            for i = 1 : 4
                J(k, i, j + 1, l) = J(k, i, 1, l) / (1 + sdtheta(j) ^ 2 * J(k, i, 1, l));
            end
        end
    end
end

fig = Figure(7, 'size', [90 90]);
styles = {'-', 'o', '+', '-'};
colors = [0 0.4 1; 0 0.4 1; 0 0.4 1; 1 0.2 0.2];
skip = [1 4 4 1 ];
start = [1 3 1 1];
set(fig, 'DefaultAxesColorOrder', repmat(colors(1 : 4, :), 2, 1))
Jmax = [5000 800];
nmax = 500;


% [A] Fisher information in the absence of input noise
subplot(2, 2, 1)
hold on
for i = 1 : 4
    plot(n(start(i) : skip(i) : end), J(start(i) : skip(i) : end, i, 1, 2), styles{i}, 'Color', colors(i, :))
end
for i = [1 4]
    plot(n, J(:, i, 1, 1), '--', 'Color', colors(i, :))
end
Jthresh = kappa ^ 2 / (sdpsi ^ 2 * mbeta ^ 2);
plot(xlim, [1 1] * Jthresh, '--k')
text(nmax, Jthresh, sprintf('%.2g', dp ./ sqrt(Jthresh) / pi * 180), 'FontSize', 14)
set(gca, 'xlim', [0 nmax], 'xtick', 0 : 100 : 500, 'ylim', [0 Jmax(1)])
xlabel('Number of neurons')
ylabel('Fisher information')
legend({'No fluctuations', 'Spatial gain', 'Feature gain', 'Attended feature'})
axis square


% [B+C] Fisher information with low+high input noise
for i = 1 : nt
    subplot(2, 2, 2 + i)
    hold on
    plot(n, J(:, 1, 1, 2), ':k')
    for j = 1 : 4
        plot(n(start(j) : skip(j) : end), J(start(j) : skip(j) : end, j, i + 1, 2), styles{j}, 'Color', colors(j, :))
    end
    for j = [1 4]
        plot(n, J(:, j, i + 1, 1), '--', 'Color', colors(j, :))
    end
    for j = 1 : 2
        Jthresh = kappa ^ 2 / (kappa ^ 2 * sdtheta(i) ^ 2 + (j - 1) * sdpsi ^ 2 * mbeta ^ 2);
        plot(xlim, [1 1] * Jthresh, '--k')
        text(nmax, Jthresh, sprintf('%.2g', dp ./ sqrt(Jthresh) / pi * 180), 'FontSize', 14)
    end
    set(gca, 'xlim', [0 nmax], 'xtick', 0 : 100 : 500, 'ylim', [0, Jmax(i)])
    xlabel('Number of neurons')
    axis square
end

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
