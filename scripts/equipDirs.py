import os

pVals = [0.75, 0.80, 0.85, 0.95]
s0Vals = [5.2,5.8]
numRuns = 30

def main():
    for p in pVals:
        for s0 in s0Vals:
            for i in range(numRuns):
                testDir = "samples_{:.2f}/{}/{}_{}/".format(p,s0,s0,i)
                os.system("cd {} && python3 ../../../3DVM_executables/sample2_5.py"
                        .format(testDir))
                os.system("cd {} && python3 ../../../3DVM_executables/ECM64.py"
                        .format(testDir))

if __name__ == "__main__":
    main()
