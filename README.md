Runs a CUDA program to calculate a statistic on the distributions of two datasets

Written to run on a SLURM cluster, to run the program:

	make run

To profile the program:

	make prof

Or:

	make vprof

and transfer the profile over to a computer with NVIDIAs visual profiler.
