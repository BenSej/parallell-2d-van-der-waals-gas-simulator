# Campus Cluster (2021)

    Campus Cluster now runs SLURM.
    To submit your job to run, enter the base directory of your project and type the following:
    
    
    #### SLURM
    ```
    sbatch scripts/batch_script.slurm
    ```
    
    That will run all benchmarks.
    
## Regression Testing
    
    There is not an easy way to unit test part3 without giving away the answers. You can test that your part3 is correct by "regression testing", where you run your code and see if its output matches known correct output.
    
    To do this, go to the base directory of the project and type:
    
    #### SLURM (Campus Cluster in 2021)
    ```
    sbatch scripts/regression_testing.slurm
    ```
    
    #### PBS/torque (Campus Cluster in 2019)
    ```
    qsub scripts/regression_testing.bash
    ```
    
    #### VMFarm, Docker, native bash shell, etc.
    ```
    bash scripts/regression_testing.bash
    ````
    
    ### Results
    This will create a file named `writeup/regression.txt`.
    The final line of this file will tell you if you pass or fail the regression test.
