# Cheet Sheet for using surface and volumetric neuroimaging data

This includes broad descriptions of the most commonly used file types Nifti, Cifti and Gifti, usage of some python libraries (mostly nibabel, nilearn), the connectome workbench (mostly for data from the human connectome project (HCP)), FSL and rarely freesurfer and AFNI.

## Nifti

Niftis are files that usually contain 3-dimensional volumetric data (except for when they are misused for surface data storage or when used as cifti files). Standard brain scan data. They also contain an affine matrix that describes how this data can be anchored into space (i.e. to align volumetric data of different granularity - or with some roation - with each other. They also have headers.

### Atalses in Volumetric MNI space

A good resource can be found on [Lead-DBS](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/cortical-atlas-parcellations-mni-space/) and ...


### Creation of a new NIFTI image

Create an empty nifti aligned to nothing:

```python
import nibabel as nib
data = np.zeros((40,40,40));
img = nib.Nifti1Image(data, np.eye(4))
```

Create a nifti aligned to an existing one:

```python
nimg = nib.load("/path/to/image.nii.gz");
newdata = ...

# with preserving header
new_image = nib.Nifti1Image(newdata, nimg.affine, header=nimg.header)

# with new header
new_image = nib.Nifti1Image(newdata, nimg.affine)
nib.to_filename("path/to/newfile.nii.gz");

```

Depending on the header size, either Nifti1 or Nifti2 images can be created (?), see [here](https://stackoverflow.com/questions/44397617/change-data-type-in-numpy-and-nibabel):
```python
hd = nimg.header
if hd['sizeof_hdr'] == 348:
  new_image = nib.Nifti1Image(newdata, nimg.affine, header=nimg.header)
else: # hdsohdr == 540
  new_image = nib.Nifti2Image(newdata, nimg.affine, header=nimg.header)
```

### Get MNI coordinate for a specific voxel

*Works only for nifti-images aligned with MNI space*: Apply the affine contained in the header:

```python
from nibabel.affines import apply_affine

nimg = nib.load("path/to/image.nii.gz");
point_xyz = [100,30,50];
# to get the value at xyz: nimg.get_fdata()[100,30,50] -> i.e. 5

mni_coords = [np.int_(np.round(x)) for x in apply_affine(nimg.affine, point_xyz)]
print("Point: ", point_xyz, "has MNI coordinates: ", mni_coords);
```

This can be reverseed using: 
```python
point_xyz = apply_affine(np.linalg.inv(nimg.affine), mni_coords)
```
but better check again.

### Get centroid for an ROI

I.e. all voxels belonging to an ROI have the value 5, then one can do:

```python

data= nimg = nib.load("path/to/image.nii.gz").get_fdata()
value=5
idx = np.array(np.where( np.round(data) == value ) ).T
centroid = [np.mean(idx[:,x]) for x in range(0,3)]

# to get useable (aka integer) voxels xyz indices, round:
centroid_round = [np.round(np.mean(idx[:,x])) for x in range(0,3)]
```


### Conversion of volumetric data (in MNI) to surface data in fsaverage space

There is multiple ways for doing so:
* there is a paper describing a good way to do it by [Wu et al., 2018 in Neuroimage](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6239990/)
* there is a [python library "regfusion"](https://pypi.org/project/regfusion/) with detailed usecase descriptions. This implements the Wu et al., 2018 mapping.
* a confusing way is given on the [freesurfer website](https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems)

This results in either (a) cifti or gifti file(s).


### HCP: "Native T1 (strutural scan)"-aligned data to surface data


```shell

# OPTIONAL:
# when the volumectric_data.nii.gz was created with FSLeyes, another step is required
# to correct some issues in the header, see
#http://yingywang.blogspot.com/2011/09/fsl-feeds-error-work-around.html
# for that 3Dcopy requires AFNI to be installed
3dcopy volumetric_data_nativeT1algined.nii.gz volumetric_data_nativeT1algined_c.nii.gz

# sometimes the axis are in the wrong orientation aparently
# so use the fsl reorientation capability:
${FSLDIR}/bin/fslreorient2std volumetric_data_nativeT1algined_c.nii.gz volumetric_data_nativeT1algined_cro.nii.gz

# Rigid Alignment with MNI (i.e. ration and resizing?)
applywarp --rel --interp=spline -i "volumetric_data_nativeT1algined_cro.nii.gz"  -r MNI152_T1_0.7mm.nii.gz --premat=/data/hcp/sub-subjectname/T1w/xfms/acpc.mat -o "volumetric_data_MNIacpc.nii.gz"

# OPTIONAL:
# if the original image was a small part of the brain (i.e. a mask or a tag)
# threshold the warped image, as the alignment introduces some kind of noise (?)
# this command requires FSL to be installed
fslmaths volumetric_data_MNIacpc.nii.gz -thr 0.6 volumetric_data_MNIacpc_clean.nii.gz

# finally do the actual volume to surface mapping
# here the target is the fsaverage_LR32k (fslr32k) space (can be sustituted as desired)
# this uses the enclosing flag; 
#-trilinear might be working better in some cases and seems to be faster

wb_command -volume-to-surface-mapping volumetric_data_MNIacpc_clean.nii.gz /data/hcp/sub-subjectname/MNINonLinear/fsaverage_LR32k/sub-subjectname.L.pial.32k_fs_LR.surf.gii surface_data.enc.32kfslr.func.gii -enclosing


```



## Cifti

Cifties can contain multiple brain structures (in contrast to Gifties?). Yet they do not contain any information about where in 3D the greyordinates (vertices) are located in, this is only contained in surf.gii files)

* .dscalar.nii(.gz), usually contains one to a hand full of scalar images, such as i.e. local cortical curvature and thickness, resulting in a shape of (2,91282) ~(n_scalars, n_vertices)
* .dtseries.nii(.gz), usually contains a time series (such as resting state) with a shape like (1200, 91282), 1200 timepoints x 91282 svoxels7vertices; 
* .dlabel.nii(.gz), contains a parcellation of shape (n_vertices) where the values correspond to a structure ID (i.e. id=4 could hypothetically correspond to the hippocampus). It should also somewhere contain the mapping from id (4) to label (hippocampus).

Internally these should look mostly similiar. For more info check out the [Layman’s guide to working with CIFTI files](https://mandymejia.com/2015/08/10/a-laymans-guide-to-working-with-cifti-files/) by Mandy Mejia. Or the [HCP FAQ on Cifti](https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ)

This format I think was created for/by the Human connectome project (HCP; and hence examples may focus on it a bit more).

### Create a new Cifti (left cortex excluding medial wall)

Each Cifti file needs an assocaited brain model (i.e. which areas are included in the cifti file, such as left hemisphere aka "LEFT_CORTEX", right hemisphere etc). This brain model we can take from an already existing file, such as the resting state run of an HCP subject:

```python
import nibabel as nib

resting_state_run = "/data/hcp/sub-100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii"
img = nib.load(resting_state_run)

# extract brain model from file
# here we just take the left hemisphere. the brain models used in the HCP are special as they exclude the medial wall
# which reduces the number of greyordinates per hemisphere from 32492 to 29696

left_cortex_excluding_medial_wall = list(img.header.matrix._mims[1].brain_models)[0].vertex_indices._indices

# create a dummy array for the greyvoxel data, and only set the vertices that are not part of the medial wall to 1
mask = np.zeros((32492)); 
np.put(mask, left_cortex_excluding_medial_wall, 1)

# create a new brain model with only the left cortex exluding medial wall
bm_left = nib.cifti2.BrainModelAxis.from_mask(mask, "LEFT_CORTEX")
```

Now we can finally create the actual image using this brain model. Multiple scalar images can be saved in a single cifti (i.e. curvature and cortical thickness, distrubution of different receptor types). Each of the scalar images needs a title:

```python
data = np.zeros((2, 29696)); # (n_scalars, n_vertices), for the current brainmodel n_vertices = 29696
# fill in / manipulate image data ...

dsnames = ["Curvature", "Thickness"];
cimg = nib.Cifti2Image(data, nib.cifti2.Cifti2Header.from_axes((nib.cifti2.ScalarAxis(dsnames), bm_leftown)))

# optionally save it to a file
cimg.to_filename("path/to/new/file.dscalar.nii(.gz)");
```

### HCP: Load resting state run Cifti and extract left cortex data

```python
resting_state_run = "/data/hcp/sub-100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii"
img = nib.load(resting_state_run)
fullbrain_tseries = img.get_fdata();
# fullbrain_tseries.shape ~ (1200, 91282), 1200 timepoints x 91282 svoxels; 

# 91282 present the whole brain, we just need the ones contained in the left cortex
# for that we can refer to the header:
structure_index = 0; # index for CORTEX_LEFT
list(img.header.matrix._mims[1].brain_models)[structure_index].brain_structure   # -> CIFTI_STRUCTURE_CORTEX_LEFT
list(img.header.matrix._mims[1].brain_models)[structure_index].index_offset    	 # -> 0 (index of first vertex contained in the structure)
list(img.header.matrix._mims[1].brain_models)[structure_index].index_count       # -> 29696 (number of greyordinates/vertices contained in the structure)
```

This resting state run is aligned by default to the FS_LR32k space. That kinda means each hemisphere is formed by ~32k greyordinates/vertices. But here we only have 29696 vertices for the left hemisphere. This is because the medial wall was excluded (as it is nessesary for the 3D mesh to exist, but does not contain any grey matter in reality). Exluding it saves storage space.

To get the resting state data for the LH only, simply do:

```python
left_ctx_tseries = fullbrain_tseries[:, 0: 0+29696];
```

The full HCP cifti greyordinate data composition looks like:

```
# LH: CORTEX_LEFT  >> N_first = 0, N_cnt = 29696 
# RH: CORTEX_RIGHT  >> N_first = 29696, N_cnt = 29716
# UL: ACCUMBENS_LEFT >> N_first = 59412, N_cnt = 135
# UR: ACCUMBENS_RIGHT >> N_first = 59547, N_cnt = 140
# ML: AMYGDALA_LEFT  >> N_first = 59687, N_cnt = 315
# MR: AMYGDALA_RIGHT  >> N_first = 60002, N_cnt = 332
# BS: BRAIN_STEM   >> N_first = 60334, N_cnt = 3472
# CL: CAUDATE_LEFT  >> N_first = 63806, N_cnt = 728
# CR: CAUDATE_RIGHT >> N_first = 64534, N_cnt = 755
# EL: CEREBELLUM_LEFT >> N_first = 65289, N_cnt = 8709  
# ER: CEREBELLUM_RIGHT  >> N_first = 73998, N_cnt = 9144
# DL: DIENCEPHALON_VENTRAL_LEFT  >> N_first = 83142, N_cnt = 706 
# DR: DIENCEPHALON_VENTRAL_RIGHT  >> N_first = 83848, N_cnt = 712
# HL: HIPPOCAMPUS_LEFT  >> N_first = 84560, N_cnt = 764
# HR: HIPPOCAMPUS_RIGHT  >> N_first = 85324, N_cnt = 795 
# PL: PALLIDUM_LEFT  >> N_first = 86119, N_cnt = 297
# PR: PALLIDUM_RIGHT >> N_first = 86416, N_cnt = 260
# AL: PUTAMEN_LEFT  >> N_first = 86676, N_cnt = 1060
# AR: PUTAMEN_RIGHT  >> N_first = 87736, N_cnt = 1010
# TL: THALAMUS_LEFT  >> N_first = 88746, N_cnt = 1288
# TR: THALAMUS_RIGHT >> N_first = 90034, N_cnt = 1248
# full : all of them >> N_first = 0, N_cnt = 91282
```

This information was taken from a [script](https://github.com/NeuroanatomyAndConnectivity/hcp_corr) by [Şeyma Bayrak](https://github.com/sheyma).


### get outline of a surface ROI

Dilate the original image (OI) to get (D). Erode original image to get (E). The output/outline image (OUT) equals then D-E. This requires the connectome workbench cli to be installed, and access to the sphere-ical surf.gii (containing postion and size of vertices in 3D) in the same surface space as the data (here: freesurfer_LR32k / fslr32k / 32k_fs_LR). Gii/Gifti files are described in the next big section.

```shell

# usage: wb_command -cifti-dilate <cifti-in> COLUMN <surface-distance> <volume-distance> <cifti-out> -left-surface /data/pt_02189/MARHCP/BrocaConn/atlases/HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.sphere.32k_fs_LR.surf.gii"

wb_command -cifti-dilate tmp921.LH29k.dscalar.nii COLUMN 2 2 tmp921.LH29k_dilated.dscalar.nii -left-surface /data/pt_02189/MARHCP/BrocaConn/atlases/HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.sphere.32k_fs_LR.surf.gii
wb_command -cifti-erode tmp921.LH29k.dscalar.nii COLUMN 2 2 tmp921.LH29k_eroded.dscalar.nii -left-surface /data/pt_02189/MARHCP/BrocaConn/atlases/HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.sphere.32k_fs_LR.surf.gii

wb_command -cifti-math "D-E" tmp921.LH29k_outline.dscalar.nii -var D tmp921.LH29k_dilated.dscalar.nii -var E tmp921.LH29k_eroded.dscalar.nii

```

### Separate Cifti files into gifti:

```shell
wb_command -cifti-separate Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_LEFT Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii

wb_command -cifti-separate Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_RIGHT Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.label.gii
```

This uses the [connectome workbench cli](). Example taken from Kathryn Mills Figshare Post [HCP-MMP1.0 projected on fsaverage](https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446).





## Gifti


* `.surf.gii` files provide the 3D brain mesh/object backbone. Contains the mapping from greyordinate vertex to its location in 3D.
* `shape.gii` or `func.gii` contain actual data, such as cortical thickness etc. (similiar to what can be contained in cifti files?)

For more info, check the posts by [Emma Robinson on Gifti/HCP surface files](https://emmarobinson01.com/2016/02/10/unofficial-guide-to-the-hcp-surface-file-formats/) and both posts by Joset (Jo) A. Etzel on [NIfTI, CIFTI, GIFTI in the HCP and Workbench](http://mvpa.blogspot.com/2014/03/nifti-cifti-gifti-in-hcp-and-workbench.html) and on [Conversion from Volumetric to Surface](https://mvpa.blogspot.com/2018/02/connectome-workbench-making-surface.html).


### Convert annoations between surface data types


*From cifti to gifti (within FS_LR32k)* 
```shell
wb_command -cifti-separate Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_LEFT Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii
```

*From FS_LR32k gifti to fsaverage164 gifti*
```shell
wb_command -label-resample Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii L.sphere.32k_fs_LR.surf.gii fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii BARYCENTRIC left.fsaverage164.label.gii
```

*From fsaverage164 gifti to fsaverage(164?) .annot-file*
```shell
mris_convert --annot left.fsaverage164.label.gii fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii lh.HCP-MMP1.annot
```

Last three examples taken from Kathryn Mills Figshare Post [HCP-MMP1.0 projected on fsaverage](https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446).

*From fsaverage .label to fsaverage gifti file*
```shell
# uses Freesurfer mris_convert
mris_convert fsaverage6-rh.label fsaverage6-rh.label.gii 
```

### Get an ROI from a Gifti parcellation

In the HCP data a parcellation is i.e. given for each subject under:
`hcp-subjectname/MNINonLinear/fsaverage_LR32k/subjectname.L.aparc.32k_fs_LR.label.gii`
(here the one in the FS_LR32k space is used)

```python
import nibabel as nib

nib.load(r"/data/hcp/hcp-subjectname/MNINonLinear/fsaverage_LR32k/subjectname.L.aparc.32k_fs_LR.label.gii")
AnatLabelsData= AnatLabels.darrays[0].data
# -> array([10, 29, 24, ..., 15, 15, 15] of shape (32492,)

# the mapping of those values to labels is given in 
AnatLabels.get_labeltable().labels

# if we were interested in the IFG pars triangularis, we would use the vertices with value 20:
AnatLabels.get_labeltable().labels[20].label 
#-> u'L_parstriangularis'

triangularis_mask = AnatLabelsData == 20;

```




### Misc: working with surface data in a 2D plane


```python


```






# HCP Processing Pipeline:

Usually human connectome data seems to be aligned to the 32k_FS_LR atlas space (but also alignments to fsaverage are available) based on the Conte69 "mesh".

Alignment as part of the HCP Pipeline works as follows: 

* PreFreesurfer Part (Volume Spaces):
[1] subject’s undistorted native volume space (rigidly “aligned” to the axes of MNI space
[2]  Standard MNI space, which is more accurately aligned by nonlinear volume registration; especially better for subcortical areas 

* FreeSurfer Part (fsaverage space):
[3] extracts white and pial cortical surfaces and transformation into parcels registered using Freesurfers standard folding-based surface registration to their surface atlas (fsaverage)

* PostFreeSurfer Part (164k and 32k_FS_LR space)
[4] Nifti(cifti?)+Gifti files for workbench; along with applying the surface
registration (to the Conte69 surface template (Van Essen et al., 2012b)), downsampling registered surfaces for connectivity analysis, creating the final brain mask, and creating myelin maps

The Pipeline is described in: Glasser, M. F., Sotiropoulos, S. N., Wilson, J. A., Coalson, T. S., Fischl, B., Andersson, J. L., ... & Van Essen, D. C. (2013). The minimal preprocessing pipelines for the Human Connectome Project. Neuroimage, 80, 105-124.


## Conte 69
Reference: Van Essen, D.C., Glasser, M.F., Dierker, D.L., Harwell, J., Coalson, T., 2012b. Parcellations and hemispheric asymmetries of human cerebral cortex analyzed on surface-based atlases. Cerebral Cortex 22, 2241-2262.

There is two versions of this atlas, a hires one using 147 vertices; and a low res one at 32k vertices per hemisphere for everyday usage. As the original atlas cannot be really retrieved anymore from their website since a technical error occured, I downloaded it from:

https://github.com/MidnightScanClub/MSCcodebase/tree/master/Utilities/Conte69_atlas-v2.LR.32k_fs_LR.wb

This LR space supposedly was first developed/described in: Van Essen, David C. “A population-average, landmark-and surface-based (PALS) atlas of human cerebral cortex.” Neuroimage 28.3 (2005): 635-662.




# HCP Data Releases
## Human HCP S900 release (Dec 2015)
## Human HCP initial 7T release (June 2016) releases
## Human HCP S1200 Structural + fMRI Atlas (March 2017; Reference)

This Human reference dataset comprises data from the Human Connectome Project (HCP) 1200 Subjects (S1200) data release; Group-average structural and functional (task) MRI data and  Selected individual-subject anatomical maps for each of the 1096 subjects having at least 1 rfMRI run 

* cites again the Glasser2013 minimal preprocessing pipelines from before
* can be seen as an extension to the conte atlas?
https://balsa.wustl.edu/reference/show/pkXDZ

MNINonLinear/ contains cortical surfaces and other data volumetrically registered to MNI152 space (using nonlinear FNIRT) followed by surface registration to Conte69 ‘164k_fs_LR’ mesh (Van Essen et al. 2012b) (via FreeSurfer fsaverage as an intermediate). *A Connectome Workbench-readable 164k_fs_LR.wb.spec file is included for quickly reading and visualizing many of these files in Workbench in the Human HCP S1200 release*

https://www.humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf 

Human.Composite_VGD11.32k_fs_LR.dlabel.nii file is derived from a composite cortical parcellation containing 52 distinct areas accurately mapped to the fs_LR atlas surface and based on architectonic or retinotopic fMRI maps (Van Essen et al. 2012). Abbreviations associated with labels of cortical areas (FRB08, OFP03, etc.) refer to the publication that defined a particular cortical area in the composite map (cf. Table 3 of Van Essen et al. 2012). http://cercor.oxfordjournals.org/content/22/10/2241.long

The Gordon333.32k_fs_LR.dlabel.nii and Gordon333_Key.txt are derived from a cortical parcellation generated using resting state functional connectivity boundary maps to identify putative borders between cortical areas Gordon et al. 2014. https://academic.oup.com/cercor/article-lookup/doi/10.1093/cercor/bhu239

RSN-networks.32k_fs_LR.dlabel.nii file displays resting state network cortical parcellation maps from Yeo et al. 2011 (7 and 17 network maps) and the Resting State network consensus communities (with and without gaps in the communities filled) from Power et al. 2011.


# Online Resource: BALSA (+ visualization assistence files ...)

Reference: DC Van Essen, J Smith, MF Glasser, J Elam, CJ Donahue, DL Dierker, EK Reid, TS Coalson, J Harwell (2016) The Brain Analysis Library of Spatial maps and Atlases (BALSA) Database.  NeuroImage (2016) PMID: 27074495; https://balsa.wustl.edu/study/show/WG33
 
 BALSA Reference is a curated repository of reference data accurately mapped to brain atlas surfaces and volumes, including various types of anatomically and functionally derived spatial maps as well as brain connectivity

Includes files for: 

Four human cortical parcellations on HCP S900 surface (left hemisphere)
Four human cortical parcellations on HCP S900 surface (montage views)
VDG11b 52-surface-mapped cortical areas (annotated)
Brodmann (1909) areas (annotated)
Gordon et al. (2016) areas (annotated)
Yeo et al. (2011) 17-network RSNs (annotated)


# HCP Atlas 

## HCP Multimodal Parcellation v1.0 (MMP1.0), Glasser et al. (Nature, 2016) 

Just describes how the parcellation is done, not which dataset?
>> Using multi-modal magnetic resonance images from the Human Connectome Project (HCP) and an objective semi-automated neuroanatomical approach, we delineated 180 areas per hemisphere bounded by sharp changes in cortical architecture, function, connectivity, and/or topography in a precisely aligned group average of 210 healthy young adults. https://pubmed.ncbi.nlm.nih.gov/27437579/

Usual File Names:
Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii



# Useful Links:

## As of October 2020

NIfTI, CIFTI, GIFTI in the HCP and Workbench: a primer 
http://mvpa.blogspot.com/2014/03/nifti-cifti-gifti-in-hcp-and-workbench.html

Getting started with Connectome Workbench 1.4.2 
https://mvpa.blogspot.com/2020/03/getting-started-with-connectome.html

Connectome Workbench: making a surface version of a volumeric image 
https://mvpa.blogspot.com/2018/02/connectome-workbench-making-surface.html

***Unofficial*** Guide to the HCP surface file formats
https://emmarobinson01.com/2016/02/10/unofficial-guide-to-the-hcp-surface-file-formats/
An important thing to recognise first about the HCP surface file format is that it has two versions of the atlas space: 164k_FS_LR and 32k_FS_LR. These atlases are regularly spaced and represent a left-right symmetric atlas developed  by Washu U in [3]. FS stands for FreeSurfer, and indicates the atlas is related to the FreeSurfer atlas fsaverage.


