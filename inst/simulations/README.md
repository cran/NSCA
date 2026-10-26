# Simulations installed with NSCA 0.5.2

These scripts reproduce the simulations referred to in NEWS.md and in the
help pages of NSCA 0.5.2. Each script has a corresponding `-result.txt` file
containing the output obtained when the release was prepared (R 4.3.3 on
Linux). The reported rates are Monte Carlo estimates and are subject to the
error stated in each script; they describe the simulated designs and may not
generalize to other designs.

| Script | What it examines | Run time |
|---|---|---|
| `calibration-0.5.2.R` | False-positive rate of each p-value when X and Y are unrelated | about 20 minutes on four cores |
| `power-0.5.2.R` | Power of the joint test and of the two supplementary tests | about 15 minutes on three cores |
| `one-side-0.5.2.R` | Designs in which one side holds and the other is absent | about 25 minutes on two cores |
| `mechanism-0.5.2.R` | Why the test of an absent corner may reject when the opposite corner is empty | about 30 minutes on one core |
| `nca-mirror-0.5.2.R` | The same tendency in the single-corner test used by NCA | about 20 minutes on one core |

Each script can be run directly from the installed package:

```r
source(system.file("simulations", "one-side-0.5.2.R", package = "NSCA"))
```

On Windows, the scripts run on a single core; this affects the run time but
not the results. None of the scripts writes files.
