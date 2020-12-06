function show_histogram(I, level)

[color1, color2, color3] = import_custom_colors();

figure()
h = histogram(I);
h.NumBins = double(max(max(I)));
h.FaceColor = [color1(1) color2(1) color3(1)];
values = h.Values;

gMin = double(min(min(I)));
gMax = double(max(max(I)));
n = round(length(values) / 2);
values1 = values(1:n);
values2 = values(n+1:end);
[~, whichbin] = max(values1);
gB = whichbin;
[~, whichbin] = max(values2);
gW = length(values1) + whichbin;

xlabel('g')
ylabel('p(g)')
axis tight
legendString = ['g_{Min} = ', num2str(gMin), ', g_{Max} = ', num2str(gMax), ...
    ', g_{B} = ', num2str(gB), ', g_{W} = ', num2str(gW)];
legend(legendString)

minV = double(min(min(I)));
maxV = double(max(max(I)));
level = level * 256;

set(gca,'xtick', minV : 50 : maxV, 'xlim', [minV, maxV])
%line([level level],[0 max(values)], 'color', [color1(2), color2(2), color3(2)], ...
%    'LineWidth', 2, 'LineStyle', '--')

%legendStringLevel = ['Threshold \tau = ', num2str(level)];
%legend(legendString, legendStringLevel)

% Optional
set(gca, 'YScale', 'log')

end