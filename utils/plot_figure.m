function quantile_freq = plot_figure(G, freq, centerbins, rep, prefix, suffix)

% plot_figure plots the degree distribution of the adjacency matrix G
% with the 75% credible intervals and returns the quantile frequencies
%
% INPUT
%   - G: observed binary adjacency matrix
%   - freq: frequencies for the counts of the degrees in each bin
%   - centerbins: bin centers
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%
%   OUTPUT
%   - quantile_freq: quantiles of the values in data for the cumulative probability




  quantile_freq = quantile(freq, [.025, .975]);
  plot_variance = @(x,lower,upper,color) fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color, 'EdgeColor', color);

  figure; hold on
  plot(centerbins, quantile_freq, 'color', [.8, .8, 1], 'linewidth', 2.5);
  ind = quantile_freq(1,:)>0;
  ha = plot_variance(centerbins(ind), quantile_freq(1,ind),quantile_freq(2,ind), [.8, .8, 1] );
  set(gca,'XScale','log')
  set(gca,'YScale','log')

  hb = plot_degree(G);
  set(hb, 'markersize', 10, 'marker', 'o',...
      'markeredgecolor', 'none', 'markerfacecolor', [1, .75, .75])

  legend([ha, hb], {'95% posterior predictive', 'Data'})
  legend boxoff
  xlim([.8, 1e3])
  box off
  set(gca,'XMinorTick','on','YMinorTick','on')

  if ~isempty(rep)
      savefigs(gcf,  [prefix 'degreepostpred' suffix], rep);
  end

end
