function cm = bluered(p)
% Blue-white-red colormap.

if ~nargin
    p = 100;
end
g = (1/p : 1/p : 1)';
gi = flipud(g);
o = ones(p, 1);
cm = [g g o; o gi gi; 0 0 0];
