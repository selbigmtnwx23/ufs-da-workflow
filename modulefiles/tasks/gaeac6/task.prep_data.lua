prepend_path("MODULEPATH", os.getenv("modulepath_spack_stack"))
load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))
load(pathJoin("nco", nco_ver))
load(pathJoin("prod_util", prod_util_ver))

prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.2/envs/ue-intel-2023.2.0/install/modulefiles/gcc/12.3.0")
load("ecflow/5.11.4")

prepend_path("MODULEPATH", os.getenv("modulepath_pymodule"))
load("python-ufs-land-da-wflow")

