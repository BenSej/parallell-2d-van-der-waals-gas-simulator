# Building

To build, create a new directory, anywhere. We will call it "BLD".
Change to that directory and type

```
cmake <path_to_source_directory>
make
```
in your shell.

This will create programs named
	* `BLD/bin/part1`
	* `BLD/bin/part2`

# Running the programs

## Running tests (locally)

On your local machine, VM, or Docker container *(not on the cluster)*,
you can run the tests we provide by executing the program
`mpirun -n 4 BLD/bin/run_tests`
(If you have problems with that, try `mpirun --launcher fork -n 4 BLD/bin/run_tests`)

Try the "--help" option to see various options.

## Running tests and benchmarks (on Campus Cluster)

See the file `README_CAMPUS_CLUSTER.md` . As always, you should submit the
appropriate script as a batch job. (Do not run intensive computations on your
login session. That is, do not run the `mpirun` commands above directly on
your login session.)

## Custom tests

See the documentation at [Google Test GitHub Project](https://github.com/google/googletest) .


# Visualization and File Format

	https://ovito.org/manual/usage.import.html#usage.import.formats
