# fMRI_prep_and_stats_(using_SPM12)
A collection of MATLAB scripts for preprocessing and statistical analysis of functional MRI data using the SPM12 toolbox.

A standard preprocessing and 1st level statistical analysis pipeline can be done stepwise using the following scripts:

`1. dicom_import.m and`  
`   dicom_import_job.m`  
The raw dicom files **(.dcm)** need to be imported as NIfTI **(.nii)** files in order to be used with SPM.  
This script imports files from a complicated folder directory structure and needs to be adapted to the  
specific directory structure of the user.

`2. convertTo4Dnifti.m`  
`   threeDto4D_job.m  `  
This step converts the previously imported single 3D NIfTI files into one 4D file (per experimental run).  
Doing this is optional and doesn't affect the the outcome of your analyses, but one or a few 4D files might  
be more managable than lots of 3D files per participant. 

`3. preprocessing.m    `  
`   preprocessing_job.m`  
Preprocessing without spatial normalisation to standard (MNI) space are performed, including *motion correction*,  
*slice-timing correction*, *co-registration to a structural image* (usually T1) and *segmentation of the structural  
image*. The temporal order of the first two steps can be swapped. **MATLAB**s parallelisation toolbox is used to  
speed up the procedure. 

`4. DARTEL_create_template.m    `  
`   DARTEL_create_template_job.m`  
The DARTEL algorithm normalizes functional images to MNI space. This procedure takes 2 steps. Here a group template  
of all participants (or participant groups with potentially different anatomy - young and old participants for example)  
is being created, alongside participant-specific **flow field images**. 

`5. DARTEL_normalize2mni.m    `  
`   DARTEL_normalize2mni_job.m`  
These scripts perform the second step of the DARTEL procedure - the actual normalization ("writing the functional into MNI  
space). This step includes the spatial smoothing of the normalized images. 

`6. check_images.m`  
Now at the latest, you should check whether the *co-registration* and *normalization* procedure have actually worked as they  
are supposed to. This script helps you to automatically go through your functional images and compare them to your structural  
image or to each other. 

`7. skull_stripping.m`  
Occasionally, something goes wrong for a certain data set when trying to co-register the mean-functional to the structural  
image, which in turn will affect the normalization step. When proper co-registration fails, one can strip the skull of the  
structual image and re-do the co-registration step. This script performs said procedure.

`8. get_headmotion_info.m`  
`   motion_censoring.m   `  
`   volt_exp_multsub.m   `  
This step is optional as well. These scripts perform various functions regarding the plotting of head motion, calculating  
the so-called Volterra expansion of motion parameters and calculating framewise discplacement values which can be used to  
flag and/or censor single images from subsequent statistical analyses. 

`9. create stim_onset_files.m`  
Before one can perform (standard univariate) statistical analyses of an experiment, one needs to create stimulus onset files  
which contain informationt to set up the **regressors of interest** of a **general linear model (GLM)**. This script performs  
this step, given that task stimulus information has already been calculated. Naturally, every experiment is different, hence  
this script needs to be customised to your own purposes, but the code provides a scaffold, for many onset file calculation  
scenarios.

`10. stats.m    `  
`    stats_job.m`  
Finally, these scripts perfom **1st level statistical analyses**. This includes, setting up a **GLM**, estimating and then  
calculating the desired statistical contrasts. The latter need - of course - be customised to your needs. Again, **MATLAB**s  
parallelisation toolbox is used to speed up the procedure.

`11. findHolesIn3DconImgs.m`  
Yet another optional step. It might happen that - depending on the **implicit masking threshold** (per default set to 0.8,  
which means 80% of the global mean of your datasets voxel intensity) - the contrast images resulting from the stats  
analysis will contain *holes* meaning parts of the brain where voxels are set to NaN. This script helps you to detect this.  
If that happens you might solve the issue by *lowering your implicit masking threshold* to something like 0.4, for example.

