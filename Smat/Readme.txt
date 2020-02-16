'../Scalc' source files in ./Smat/

'fopen' opens the following files.
line	read/write	function			file name
45		r			main				argv[1]	= ./input.txt
55		r			main				'f_prefix'_Med'p_target'.dat
121		w			Out_Smat			'f_prefix'_Smatrix.dat
524		r			set_Phi				Xi_file
558		r			set_tilPhi			Xi_file
608		r			get_RealNum05omega	file_name = 'f_prefix'_bMode.dat 'f_prefix'_tMode.dat
844		r			set_V				data_name

main:
	235 data_name = 'f_prefix'_MedBj.dat 'f_prefix'_Medj.dat 'f_prefix'_MedTj.dat

initial_calc:
	190	Xi_file = 'f_prefix'_bXi.dat
	191	Xi_file = 'f_prefix'_tXi.dat

input_file0:
	666	'f_prefix' = ./Data/joint
	671	'p_target' = B7
	673	'p_target' = T1
