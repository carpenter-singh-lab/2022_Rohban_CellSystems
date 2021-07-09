figure(99)
clf

plot(newtab.Broad(newtab.isrndm == 1 & newtab.medTLR5 > 1), log(newtab.AdjP(newtab.isrndm == 1 & newtab.medTLR5 > 1)), 'o', 'markersize', 9, 'color', [0.7 0.7 0.7], 'markerfacecolor', [0.7 0.7 0.7]); hold on;
plot(newtab.Broad(newtab.isrndm == 1 & newtab.medTLR5 < 1), log(newtab.AdjP(newtab.isrndm == 1 & newtab.medTLR5 < 1)), 'o', 'markersize', 9, 'color', [0.3 0.3 0.3], 'markerfacecolor', [0.3 0.3 0.3]);
plot(newtab.Broad(newtab.isrndm == 0 & newtab.medTLR5 > 1), log(newtab.AdjP(newtab.isrndm == 0 & newtab.medTLR5 > 1)), 'o', 'markersize', 9, 'color', [1 0.5 0.5], 'markerfacecolor', [1 0.5 0.5]); hold on;
plot(newtab.Broad(newtab.isrndm == 0 & newtab.medTLR5 < 1), log(newtab.AdjP(newtab.isrndm == 0 & newtab.medTLR5 < 1)), 'o', 'markersize', 9, 'color', [0.6 0 0], 'markerfacecolor', [0.6 0 0]);

line(get(gca, 'xlim'), [log(q(end)) log(q(end))], 'linestyle', '--', 'color', [0.5 0.5 0.5])
xlabel('Cell Painting Correlation'); ylabel('FDR adjusted p-value (log)'); axis([- 0.6 0.7 - 41 - 5])
set(gca, 'ytick', - 40:10:0, 'xtick', - 0.5:0.5:0.5, 'fontsize', 13)