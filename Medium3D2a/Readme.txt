'../Med3D2' source files in ./Medium3D2a/

"Med3D2a" is an updated version of "Med3D".
We can set smoothing factors for epsilon and mu on lattice points.

The following is previous Readme for "Med3D".

'fopen' opens the following files.
line	function		file name
33 r	main			argv[1]	= input.txt
54 r	main			MediumTable
195 w	output_file		output_name = %s_Med%%.dat
295 r	get_file_info		input_name
326 r	input_data_file		input_name
413 r	input_data_file0	input_nameB = %sB0.dat
429 r	input_data_file0	input_name = %s0.dat
446 r	input_data_file0	input_nameT = %sT0.dat

input_data_file:
		input_name = %s%%.dat

input_file1:
		MediumTable = ./trans_med.txt

