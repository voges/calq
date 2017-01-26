function generate_P_R_F_plots(avg_data, avg_size, file_name)

    metric = {'Recall', 'Precision'};
    h = figure;
    for j = 1:length(metric)
        P = avg_data(:,j:3:end);
        plot_joint_PRF(avg_size, P, metric{j}, j)

    end

    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    print(h, '-dpdf', [file_name '_Recall_Precision_versus_bits_per_QV.pdf']);

end

