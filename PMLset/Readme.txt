'../PMLgen' source files in ./PMLset/

'fopen' opens the following files.
line	read/write	function			file name
39		r			main				argv[1] = input.txt
127		r			input_u2			'Yee'_Medn.dat
168		r			output_PML			'Yee'_Medn.dat 
169		w			output_PML			'Yee'_PMLn.dat renamed to 'Yee'_Medn.dat
362		r			input_u				u0u1file
397		r			input_PMLn			inputname = 'Yee'_Med0.dat
412		r			input_PMLn			PMLfile
529		r			input_L				inputname = 'Yee'_Med0.dat
611		r			input_N_str			PMLfile

input_filename:
	555 u0u1file = ./non_uniform.dat
	558	'Yee' = ./Data/joint
	561 PMLfile = ./PML_str.txt

