function [data_opt_t, t_opt]  = compute_optimal_t(data, DiffOp, varargin)

t_max = 32;
make_plot = true;
th = 1e-3;
data_opt_t = [];

if ~isempty(varargin)
    for j = 1:length(varargin)
        if strcmp(varargin{j}, 't_max')
            t_max = varargin{j+1};
        end
        if strcmp(varargin{j}, 'make_plot')
            make_plot = varargin{j+1};
        end
        if strcmp(varargin{j}, 'th')
            th = varargin{j+1};
        end
    end
end

data_prev = data;
if make_plot
    error_vec = nan(t_max,1);
    for I=1:t_max
        disp(['t = ' num2str(I)]);
        data_curr = DiffOp * data_prev;
        error_vec(I) = procrustes(data_prev, data_curr);
        if error_vec(I) < th && isempty(data_opt_t)
            data_opt_t = data_curr;
        end
        data_prev = data_curr;
    end
    t_opt = find(error_vec < th, 1, 'first');
    
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
else
    for I=1:t_max
        disp(['t = ' num2str(I)]);
        data_curr = DiffOp * data_prev;
        error = procrustes(data_prev, data_curr);
        if error < th
            t_opt = I;
            data_opt_t = data_curr;
            break
        end
        data_prev = data_curr;
    end
end

disp(['optimal t = ' num2str(t_opt)]);









