clear all;
run 'data_container.m'
fields = fieldnames(s)';

for i = 1: length(fields)
    generate_P_R_F_plots_indv(s.(fields{i}).diff_data,s.(fields{i}).size);
end

