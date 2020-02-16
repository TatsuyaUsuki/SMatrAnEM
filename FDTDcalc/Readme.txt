'../FDTD' source files in ./FDTDcalc

'fopen' opens the following files.
line	read/write	function	file name
50		r			main		argv[1]	= ./input.txt
560		r			set_Xi		Xi_file		
601		r			set_tilXi	Xi_file
635		r			get_RealNum	file_name
668		r			input_data	data_name
729		r			input_Lf	Yee_Med0.dat
743		r			input_Lf	Yee_bMode.dat

main:
	59	Yee = f_prefix

S_calc:
	434 file_name = mode_file

set_incident:
	483	Xi_file = f_prefix_bXi.dat
	484 mode_file = f_prefix_bMode.dat
	486	Xi_file = f_prefix_tXi.dat
	487 mode_file = f_prefix_tMode.dat

matdata_file:
	716	data_name = f_prefix_MedBj.dat 
	717	data_name = f_prefix_Medj.dat
	718	data_name = f_prefix_MedTj.dat

input_file0:
	765 f_prefix = ./Data/joint

