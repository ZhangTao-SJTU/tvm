# A python script to automate a sample run of the code. 
# 
# Executing this script should 
# 1. Recompile the code (or compile it for the first time)
# 2. Create a new directory called "sample_run/"
# 3. Copy sample_conf to sample_run/conf
# 4. Execute scripts/tvm/sample2_5.py and scripts/tvm/ECM64.py to create the initializations sample_run/sample.topo and sample_run/ECM.topo
# 5. Run the code given the conf and initializations in sample_run/
import os
if not os.path.exists("build/"):
    os.system("chmod +x compile.sh && ./compile.sh")
if os.path.exists("sample_run/"):
    os.system("rm -rf sample_run/")
os.mkdir("sample_run/")
os.system("cp sample_conf sample_run/conf")
os.system("cd sample_run/ && python3 ../scripts/tvm/sample2_5.py sample_run")
os.system("cd sample_run/ && python3 ../scripts/tvm/ECM64.py sample_run")
os.system("cd sample_run/ && ../build/tvm")