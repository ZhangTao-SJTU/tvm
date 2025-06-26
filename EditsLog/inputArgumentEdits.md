1. Use input argument to make the code run only a single stage of FIRE minimization

After the initializations in tvm.cpp, add the following block:
if (argc > 1 && string(argv[1]) == "FIRE_only") {}