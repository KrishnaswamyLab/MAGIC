function [t_opt, error_vec] = compute_optimal_t_fast(data, DiffOp, varargin)

t_max = 20;
make_plots = true;
th = 0.001;

if ~isempty(varargin)
    for j = 1:length(varargin)
        if strcmp(varargin{j}, 't_max')
            t_max = varargin{j+1};
        end
        if strcmp(varargin{j}, 'make_plots')
            make_plots = varargin{j+1};
        end
    end
end

error_vec = nan(t_max,1);
data_prev = data;
for I=1:t_max
    data_curr = DiffOp * data_prev;
    error_vec(I) = procrustes(data_prev, data_curr);   
    data_prev = data_curr;
end

t_opt = find(error_vec < th, 1, 'first');

disp(['optimal t = ' num2str(t_opt)]);

if make_plots
    figure;
    hold all;
    plot(1:t_max, error_vec, '*-');
    plot(t_opt, error_vec(t_opt), 'or', 'markersize', 10);
    xlabel 't'
    ylabel 'error'
    axis tight
    ylim([0 ceil(max(error_vec)*10)/10]);
    plot(xlim, [th th], '--k');
    legend({'y' 'optimal t' ['y=' num2str(th)]});
    set(gca,'xtick',1:t_max);
    set(gca,'ytick',0:0.1:1);
end


