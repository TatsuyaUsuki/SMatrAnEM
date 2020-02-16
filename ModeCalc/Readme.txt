'../Mcalc' source files in ./ModeCalc/

'fopen' opens the following files.
line	read/write	function	file name
24	r		main		argv[1] = ./input.txt
133	r		main_calc	data_name
148	w		main_calc	Xi_name
161	w		main_calc	out_name


main:
	36	data_name = 'f_prefix'_Med'p_cross'.dat
	40	out_name = 'f_prefix'_bMode.dat
	42	Xi_name	= 'f_prefix'_bXi.dat
	46	out_name = 'f_prefix'_tMode.dat
	48	Xi_name	= 'f_prefix'_tXi.dat

input_file0:
	67	'f_prefix' = ./Data/joint
	72	'p_cross[0]' =  B7
	75	'p_cross[BUFSIZE]' = T1
