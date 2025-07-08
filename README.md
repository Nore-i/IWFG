You need to run `stk_test.m` in the root directory once to compile them into mex files compatible with your operating system, and then in ./misc/pareto, there is an IWFG_test.m, and it demonstrates how to use these two algorithms.
`stk_dominatedhv` returns the hypervolume of a set. The underlying algorithm is "WFG". As I said, previously you have to call it twice to compute the exclusive hypervolume of a single point.
`stk_IWFG` returns the bottom-k hypervolume contributor of a set. The underlying algorithm is "IWFG".
