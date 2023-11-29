# online_study

# IAC_tDDM_CPL.R

## Modification from example file

1. RT in our study is in ms, was transformed to seconds
2. fit function (fitSub) now generates a csv every time a single subject fit is completed, the label is context + subject ID

| **var** |                         **info**                         | Lower boundary | Upper boundary |
|:-------:|:--------------------------------------------------------:|----------------|----------------|
|   d_v   | weighting for attribute 1 (taste)                        |       -2       |        2       |
|   d_h   | weighting for attribute 2 (health)                       |       -2       |        2       |
|  thres  | threshold                                                |       0.6      |        3       |
|   nDT   | non-decision time                                        |      0.01      |        1       |
|   tHin  | time Health in (relative start time of health vs. taste) |       -1       |        1       |
|   bias  | starting point bias                                      |       -1       |        1       |
