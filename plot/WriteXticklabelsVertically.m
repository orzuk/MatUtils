function WriteXticklabelsVertically(S, N)
% Writes X tick labels S on axes vertically (from Tal Shay).
% N is the number of rows in the current axes
nS = length(S);
set(gca, 'xtick', 1:nS, 'xticklabel', []);
for i = 1:nS
    text(i,N,S{i}, 'rotation', 270, 'interpreter', 'none');
end
