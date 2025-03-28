The file for the computation of the largest Lyapunov exponents, winding numbers and average infectious can be executed by running

```python3 SIRS-s-M_LLE_W_avI_lnPatches.py $OUTPUT_DIR $PARAMETER_DIR $TASKID_ZEROBASED $PROGRESSBOOL```

, where the output will be written in $OUTPUT_DIR, a 'Parameter.yml' file is required in the directory $PARAMETER_DIR, and the task id  $TASKID_ZEROBASED needs to range over all patch numbers (an integer 0,1,...,N-1 where N is the number of patches). By specifying $PROGRESSBOOL (true/false) progress-feedback will be returned or not.
After all patch numbers are computed, the full dataset is combined by running
```python3 CombinePatches.py  $OUTPUT_DIR $PARAMETER_DIR```