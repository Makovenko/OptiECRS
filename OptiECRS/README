OptiECRS is a routine to perform optimization of Extended Cauchy Reed-Solomon codes.
Supported methods:
	BB:    Branch And Bound
	RH:    Random Heuristic
	IP:    Mixed Integer Programming
	IPAlt: Mixed Integer Programming (Alternative Formulation)

=== BUILD ===
You will need Gurobi solver library and headers to use be able to build this code.
In Makefile, adjust GUROBIINC, GUROBILIB and GUROBICAL to point to correct Gurobi headers
and shared objects.

After that, the code should compile with:
./make

You will need C++17 compliant compiler.

=== RUN ===
To run the code, use ./OptiECRS.

Follwing is a sample interaction:

$: ./OptiECRS 
	Galois field w: 4
	Number of rows (informational packets): 6
	Number of cols (parity packets): 2
Initial:
Extended Cauchy Matrix (4, 6, 2):
	   7   6
	   6   7
	  13  11
	  11  13
	   9  14
	  14   9
Number of ones: 112
Extended Cauchy Matrix (4, 6, 2) generators:
	 x: (   6,    7)
	 r: (   1,    1)
	 y: (   0,    1,    2,    3,    4,    5)
	 s: (   1,    1,    1,    1,    1,    1)

	Timelimit (seconds): 60
	Algorithm (BB, RH, IP or IPAlt): BB
OptimizeBranchAndBound: found a new solution with cost 57
Extended Cauchy Matrix (4, 6, 2):
	   4   9
	   1   4
	   4   1
	   1   9
	   9   1
	   1   1
Number of ones: 57
Extended Cauchy Matrix (4, 6, 2) generators:
	 x: (   0,    1)
	 r: (   1,    3)
	 y: (  12,   11,    4,    6,    2,    9)
	 s: (   5,   11,    3,    6,    1,    9)

Time taken by algorithm: 0.00031085 seconds.
Extended Cauchy Matrix (4, 6, 2):
	   4   9
	   1   4
	   4   1
	   1   9
	   9   1
	   1   1
Number of ones: 57
Extended Cauchy Matrix (4, 6, 2) generators:
	 x: (   0,    1)
	 r: (   1,    3)
	 y: (  12,   11,    4,    6,    2,    9)
	 s: (   5,   11,    3,    6,    1,    9)
