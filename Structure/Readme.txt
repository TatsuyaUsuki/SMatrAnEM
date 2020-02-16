'../StrMap' source files in ./Structure/

'fopen' opens the following files.
line	function		file name
37 r	main			argv[1]	= ./input.txt
58 r	main			xi_file
95 r	main			str_file
213 r	input_file1		xi_file
602 w	output_file		buf = Yee%d.dat YeeB%d.dat YeeT%d.dat


input_file1:
	198	xi_file = ./non_uniform.dat
	203	str_file = ./3d_str.txt
	205	Yee = Prefix = ./Data/joint

