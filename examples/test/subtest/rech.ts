READ_LUMPREM_OUTPUT_FILE lr_def.out 2
#  my_name     LUMPREM_name      divide_by_delta_t?

	rech		net_recharge		div_delta_t
	pumping		gw_withdrawal		div_delta_t


READ_LUMPREM_OUTPUT_FILE lr_ghj.out 2
#  my_name     LUMPREM_name      divide_by_delta_t?

	rech		net_recharge		div_delta_t
	pumping		gw_withdrawal		div_delta_t


# Begin writting MF6 time series

#	ts_name		scale		offset		mf6method

	rech		1		0		linearend
	pumping		1		0		linearend
