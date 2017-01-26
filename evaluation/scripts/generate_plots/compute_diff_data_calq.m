function diff = compute_diff_data_calq(input_data)

diff = input_data - repmat(input_data(1,:),size(input_data, 1),1) ;

end

