function generate_P_R_F_plots_indv(data, bpq)

filter_value={'90','99','99.9','100'};
metric = {'Recall', 'Precision'};

for j=1:length(filter_value)
    figure;
    for i = 1:length(metric)
        P = data(:,i:3:end);
        plot_indv_P(bpq(2:end,:), P(:,j), filter_value{j}, metric{i}, i);
    end
end

