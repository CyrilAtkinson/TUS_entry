# Optimal TUS transducer position

This code aims to determine the best position of a neuromodulation ultrasound transducer over the head.

It requires 2 nifti files as input: the scalp and the target, both as binary maps

## Platform

This code was developed and tested on Windows 11 (intel i7; 16Go RAM) using R (version 4.1.1) 

## Dependencies
```
[cluster] ([https://cran.r-project.org/web/packages/cluster/index.html])
[cowplot] ([https://cran.r-project.org/web/packages/cowplot/index.html])
[lattice] ([https://cran.r-project.org/web/packages/lattice/index.html])
[oro.nifti] ([https://cran.r-project.org/web/packages/oro.nifti/index.html])
```

## Instructions

Install dependencies (see above). 

Copy and paste the TUS_entry R function

Load the mandatory files as follows:
```
target = readNIfTI("PATH_TO_YOUR_DATA", reorient=F)
scalp = readNIfTI("PATH_TO_YOUR_DATA", reorient=F)
```

If needed, load the optional files as follows:
```
t1 = readNIfTI("PATH_TO_YOUR_DATA", reorient=F)  # (default:NULL)
exclu = readNIfTI("PATH_TO_YOUR_DATA", reorient=F) # (default:NULL)
```

If needed, change the following parameters:
```
output = "TEXT" # (default:NULL)
visual_confirm = T/F # (default:F)
resolution = NUMERICAL_VALUE # (default:0.33cm)
voxel_size = NUMERICAL_VALUE # (default:0.1cm)
minimal_distance = NUMERICAL_VALUE # (default:4.24cm)
maximal_distance = NUMERICAL_VALUE # (default:7.46cm)
transducer_size = NUMERICAL_VALUE # (default:6.4cm)
kplan_offset = NUMERICAL_VALUE # (default:1.082cm)
```

Run the function as follows:
```
OBJECT_NAME = TUS_entry(target=target, scalp=scalp) # add others parameters if changes are needed
```

Save required outputs as follows:
```
writeNIfTI(OBJECT_NAME$MRI_Final_neuronav, "PATH_TO_YOUR_DATA", gzipped=F) # to save the nifti file for neuronavigation software
writeNIfTI(OBJECT_NAME$MRI_Final_validation, "PATH_TO_YOUR_DATA", gzipped=F) # to save the transducer as a binary map for 3D validation
sink("PATH_TO_YOUR_DATA) ; cat(OBJECT_NAME$report) ; sink() # to save the report
```

## Note

The input files should all have the same orientation, and dimensions and be aligned.
The outputs will be provided in the same orientation as the input files.
The scalp mask could be obtained with AFNI as follows:
```
3dAutomask -clfrac 0.4 -prefix OUTPUT INPUT # the clfrac could be modified and validate visually
```

Mango ([https://mangoviewer.com/mango.html]) could be usefull to generate the ROI targets and for the 3D visual validation. 

## Troubleshooting

cyril.atkinsonclement@gmail.com

## Citing this work

The full explanation and purpose of the code is described here: https://doi.org/10.1016/j.neurom.2024.06.496

If you use this code, please acknowledge us by citing the above paper and the repository.


Feedback welcome at cyril.atkinsonclement@gmail.com
