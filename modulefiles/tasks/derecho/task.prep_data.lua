prepend_path("MODULEPATH", os.getenv("modulepath_spack_stack"))
load(pathJoin("stack-oneapi", stack_intel_ver))
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

load(pathJoin("nco", nco_ver))
load(pathJoin("prod_util", prod_util_ver))

load("mkl/2024.2.2")
load("ecflow/5.11.4")

prepend_path("MODULEPATH", os.getenv("modulepath_pymodule"))
load("python-ufs-land-da-wflow")

