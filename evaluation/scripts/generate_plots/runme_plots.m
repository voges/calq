clear all;
run 'data_container.m'
fields = fieldnames(s)';

avg_data_HS01 = (s.(fields{1}).diff_data + s.(fields{2}).diff_data)/2;
avg_size_HS01 = (s.(fields{1}).size + s.(fields{2}).size)/2;

avg_data_garvan = (s.(fields{3}).diff_data + s.(fields{4}).diff_data)/2;
avg_size_garvan = (s.(fields{3}).size + s.(fields{4}).size)/2;

avg_data_iontorrent = (s.(fields{5}).diff_data + s.(fields{6}).diff_data)/2;
avg_size_iontorrent = (s.(fields{5}).size + s.(fields{6}).size)/2;

generate_P_R_F_plots(avg_data_garvan, avg_size_garvan(2:end,:), 'H12')
generate_P_R_F_plots(avg_data_HS01, avg_size_HS01(2:end,:), 'H01')
generate_P_R_F_plots(avg_data_iontorrent, avg_size_iontorrent(2:end,:), 'H11')

