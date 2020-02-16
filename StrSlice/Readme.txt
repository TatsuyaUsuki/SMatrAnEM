This 'Readme.txt' explains '../YeeSlice' version 18.09.17.
Option '-v' or '--version' shows the version and other information.

The '../YeeSlice' can create 2D cross-section data from 3D medium data to check the 3D settings.

Source file: ../StrSlice/YeeSlice.c
'fopen' in the file opens the following files.

line	read/write	function	file name
25	r		main		argv[1] = ./input.txt
161	r		input_file2	file_name = 'Yee'_MedB0.dat
170	r		input_file2	ini_file = ./input.txt after the 1st '=' in a line beginning with '# info' file_name.
211	r		Case_u2		file_name = 'Yee'_Medj.dat
214	w		Case_u2		's_name'_0.dat
215	w		Case_u2		's_name'_1.dat
216	w		Case_u2		's_name'_2.dat
280	r		get_u		non_uniform
304	r		get_u0r0	file_name is defined in Case_u0
364	w		Case_u0		's_name'_0.dat
365	w		Case_u0		's_name'_1.dat
366	w		Case_u0		's_name'_2.dat
393	r		Case_u0		file_name = 'Yee'_Medj.dat
469	w		Case_u1		's_name'_0.dat
470	w		Case_u1		's_name'_1.dat
471	w		Case_u1		's_name'_2.dat
502	r		Case_u1		file_name = 'Yee'_Medj.dat
565	r		get_u2x2	file_name is defined in Case_u0

The above string arrays are set as follows.

input_file1:
	73	's_name' = ./SliceData/Slice0, ...etc.
		's_name' are parameters after the 3rd '=' in lines with prefix 'Slice', 'slice' or 'SLICE'.
	
input_file2:
	168	non_uniform = ./non_uniform.dat
		non_uniform is a parameter after the 2nd '=' in a line beginning with '# info' in file_name.
		The original is after the 1st '=' of a line with prefix 'xi2u' in argv[1].

	174	u2r2x_name = ./urx.dat
		non_uniform is a parameter after the 1st '=' of a line with prefix 'u2r2x' in ini_file.

	175	outer_name = ./outer.dat
		outer_name is a parameter after the 2nd '=' of a line with prefix 'outer' in ini_file.

Case_u0:
	350	file_name = u2r2x_name
	352	file_name = outer_name
