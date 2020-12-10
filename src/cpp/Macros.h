#ifndef MACROS_H
#define MACROS_H

// #define OPENMP              // option to compute smooth particle distributions with OpenMP (assumes openmp support)
								// todo: use the standard macro _OPENMP instead

// #define PRINT_PARAMETERS    // option to print the runtime parameters (useful for debugging)

#define FLAGS				// option to print warnings during runtime

// #define PRINT_PROCESSOR		// option to print the processor you're using (for benchmarking)

#define ABORT_GSL				// turn on gsl abort if error in gsl matrix solver (program quits)

#define MONITOR_FAMOD			// monitor breakdown of famod and reconstruction of anisotropic variables (turn off for OpenMP or real runs)

// maybe I should make a JETSCAPE macro...

// need to write in dependencies?


#endif