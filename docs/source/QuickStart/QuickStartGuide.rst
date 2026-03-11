.. _QuickStart:

***************************
Running the UFS DA Workflow
***************************

This chapter provides instructions for building and running the Unified Forecast System (:term:`UFS`) DA Workflow System.

.. include:: ../doc-snippets/gcblizzard-desc.rst

.. attention::
   
   These steps are designed for use on :ref:`Level 1 <LevelsOfSupport>` systems (e.g., Ursa, Orion, Hercules, Gaea-C6 and Derecho.)

.. _GetCode:

Get Code
***********

Clone the Land DA workflow repository. To clone the ``develop`` branch, run:

.. code-block:: console

   git clone -b develop --recursive https://github.com/ufs-community/ufs-da-workflow

.. _BuildandRun:

Build and Run the DA Workflow
******************************

#. Navigate to the ``sorc`` directory.

   .. code-block:: console

      cd ${BASEDIR}/ufs-da-workflow/sorc

#. Run the build script ``app_build.sh``:

Note: Running and building the UFS Da Workflow options below will have different run styles based off the use of the JEDI-Bundle/GDAS and Workflow components and where '[APP]' is `S2SWA`, `S2SWAL`, `NG-GODAS`, `ATML`, or `ATM`.

#. Option 1: Workflow components: YES, JEDI-bundle: NO

   .. code-block:: console

      ./app_build.sh -a=[APP]

#. Option 2: Workflow components: YES, JEDI-bundle: YES

   .. code-block:: console

      ./app_build.sh -a=[APP] --jedi=bundle

#. Option 3: Workflow components: NO, JEDI-bundle: YES

   .. code-block:: console

      ./app_build.sh -a=[APP] --jedi=bundle-only

#. Option 4: Workflow components: YES, GDAS App: YES

   .. code-block:: console

      ./app_build.sh -a=[APP] --jedi=gdas

#. Option 5: Workflow components: NO, JEDI-bundle: YES

   .. code-block:: console

      ./app_build.sh -a=[APP] --jedi=gdas-only

.. _load-env:

Load the Workflow Environment
******************************

To load the workflow environment:

.. code-block:: console
	
   cd ..
   module use modulefiles
   module load wflow_[workflow_manager]_[machine]

where `[workflow_manager]` is `ecflow`, `rocoto`, or `none`, and `[machine]` is `gaeac6`, `hercules`, `orion`, `ursa`, or `derecho`.

Modify the Workflow Configuration YAML
**************************************

Copy the experiment settings into ``config.yaml``:

.. code-block:: console

   cd ${BASEDIR}/ufs-da-workflow/parm
   cp config_samples/config.[sample_case].yaml config.yaml

Users will need to configure the ``account`` variable in ``config.yaml``.

   * ``ACCOUNT:`` A valid account name. Most NOAA :term:`RDHPCS` systems require a valid account name; other systems may not (in which case, any value will do).
 

.. attention:: 

   When regenerating an experiment from the same or similar ``config.yaml`` file, if the ``EXP_CASE_NAME`` remains the same, the old experiment directory with that name will be renamed with the ``*_old`` suffix, and the new experiment directory will use ``EXP_CASE_NAME``. However, the ``envir`` directory will **NOT** be regenerated unless the ``envir`` parameter is given a new name. If it keeps the same name, the previous ``ptmp/<envir>`` directory and everything in it will remain (rather than being renamed), and the experiment will continue from where it left off using the files from the previous directory. This can be helpful in certain cases but detrimental in others, so users need to make a conscious choice based on their use case. 

.. _generate-wflow:

Set Up and Launch the Workflow
*******************************

Generate the experiment directory by running: 

.. code-block:: console

   ./setup_wflow_env.py 

Note: Two options for running the workflow based off the desired WORKFLOW_MANAGER ('ecflow' or 'rocoto'.)

#. Option 1: WORKFLOW_MANAGER: ecflow

.. code-block:: console

   cd ../../exp_case/ecf_server
   ./start_server.sh
   cd ../[EXP_CASE_NAME]/ecf
   ./begin_suite.sh
   ecflow_ui &

Where `[EXP_CASE_NAME]` is specified in the configuration file `config.yaml`.

#. Option 2: WORKFLOW_MANAGER: rocoto

.. code-block:: console

   cd ../../exp_case/[EXP_CASE_NAME]
   ./automate_launch_script.py -i [time interval in seconds]

Where the default value of `[time interval in seconds]` is 30. This means that the launch script `launch_rocoto_wflow.sh` is submitted every 30 seconds.


Check the result and log files
*******************************

- `com_dir`: symlink to the directory containing the result files
- `log_dir`: symlink to the directory containing the log files
- `tmp_dir`: symlink to the directory containing the working directories



