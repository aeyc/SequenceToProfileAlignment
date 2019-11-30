Ayca Begum Tascioglu
21600907

To run the code:

1) Please go to the directory of where source code(alignSeqToProfile.py) located.
 
2) Add profile file which includes multiple aligned sequences(aligned_sequences.aln), that you want to use, in the same file with the source code. Also, you must add the sequence file(seq.fasta) that you want to align to profile.

3) After, in the terminal go to the directory where alignSeqToProfile.py located and desired input files( i.e sequ.fasta and aligned_sequences.aln etc) are located.

4) type 'make'.

5) Program asks for the files which profile and sequence are located in.
	Type the name of aln file which includes multiple aligned patterns:

		>>aligned_sequences.aln

	Type the name of sequence file:

		>>seq.fasta

6) Program asks for the match, mismatch and the gap score.

	Type match score

		>>1

	Type mismatch score

		>>-1

	Type gap score

		>>-2

6) The program will run and creates seq.aln as the output file. In seq.aln, the profile(given multiple sequence alignment) and the new aligned sequence are given.

