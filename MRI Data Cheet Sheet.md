# Cheet Sheet for using surface and volumetric neuroimaging data

This includes broad descriptions of the most commonly used file types Nifti, Cifti and Gifti, usage of some python libraries (mostly nibabel, nilearn), the connectome workbench (mostly for data from the human connectome project (HCP)), FSL and rarely freesurfer and AFNI.

This makes mainly use of 
* python (specifically [Nilearn](https://nilearn.github.io/stable/index.html "Nilearn")/[Nibabel](https://nipy.org/nibabel/ "Nibabel")), 
* the [connectome workbench](https://humanconnectome.org/software/connectome-workbench "connectome workbench") CLI ([wb_command function reference](https://humanconnectome.org/software/workbench-command "wb_command function reference")) and 
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki "FSL") ([fslmaths commands](https://mandymejia.com/fsl-maths-commands/ "fslmaths commands"))

For editing of this file, use i.e. [pandao MD editor](https://pandao.github.io/editor.md/) or [stackedit](https://stackedit.io/app#)


## File Types Overview

**Standard formats**
* Nifti: .nii, .nii.gz
* Cifti: dscalar.nii, dlabel.nii, dtseries.nii, (plabel.nii)
* Gifti data: shape.gii,  func.gii, label.gii 
* Gifti mesh: surf.gii*

**Freesurfer formats**
* lh.sphere
* ?h.sulc - subject's convexity data
* ?h.smoothwm surface
* ...
* subj01_parc_masks.mgz - volumetric parcellation?
* lh.XXXX.annot - annotation/parcellation mapped to an individual (or refrence template), can also contain a color mapping, see [here](https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles#Annotation)
* XXXColorLUT.txt - label name to RGBA color mapping (weirdly enough, it doesnt map label iDs), see [here](https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles#Annotation)
* lh.XXXX.gcs - atlas file based on multiple individual annotations along with geometric (sulci/gyri) data that allows freesurfer automatically map this parcellations to individuals

See [freesurfer parcellation overview](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) and the related [automatic surface labelling process](https://surfer.nmr.mgh.harvard.edu/fswiki/SurfaceLabelAtlas)
https://github.com/binarybottle/mindboggle_sidelined/blob/master/freesurfer.py

See also the different freesurfer spaces in the end of this document. 

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

# alternative way to create new image of same class
new_image = nimg.__class__(newdata, nimg.affine, header=nimg.header)

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

### MNI aligned volumetric NIFTI -> HCP FS_LR_32k

```shell

#wb_command -volume-to-surface-mapping <input> <map_target_surf> <output> [-method]
wb_command -volume-to-surface-mapping juelich_cyto_MNI152_2009C_nA.nii.gz atlases\surf_fsLR\reference\HCP_S1200\S1200.L.pial_MSMAll.32k_fs_LR.surf.gii juelich_cyto_MNI152_2009C_nA.enc.32kfslr.L.func.gii -enclosing

wb_command -volume-to-surface-mapping juelich_cyto_MNI152_2009C_nA.nii.gz atlases\surf_fsLR\reference\HCP_S1200\S1200.R.pial_MSMAll.32k_fs_LR.surf.gii juelich_cyto_MNI152_2009C_nA.enc.32kfslr.R.func.gii -enclosing

# -enclosing does not average voxels, so can be used i.e. for transformation of atlases
# - trilinear uses interpolation; -cubic & -ribbon-constrained are even more complicated

```

### HCP: Functional Image aligned to T1 (structural scan) in subject native volumetric space -> Subject specific FS_LR_32k

"Native T1 (strutural scan)"-aligned data to surface data
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

# Rigid Alignment with MNI (i.e. rotation and resizing?)
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
* .dtseries.nii(.gz), usually contains a time series (such as resting state) with a shape like (1200, 91282), 1200 timepoints x 91282 svoxels/vertices; 
* .dlabel.nii(.gz), contains a parcellation of shape (n_vertices) where the values correspond to a structure ID (i.e. id=4 could hypothetically correspond to the hippocampus). It should also somewhere contain the mapping from id (4) to label (hippocampus).

Internally these should look mostly similiar. For more info check out the [Layman’s guide to working with CIFTI files](https://mandymejia.com/2015/08/10/a-laymans-guide-to-working-with-cifti-files/) by Mandy Mejia. Or the [HCP FAQ on Cifti](https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ)

This format I think was created for/by the Human connectome project (HCP; and hence examples may focus on it a bit more).

### Look at what is contained within a cifti

Scalar file:
```python
import nibabel as nib
cifti=nib.load(r"MNINonlinear\fsaverage_LR32k\100206.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii" ) 
# data within cifti of shape (1, 64984)
cifti.header.get_axis(0) # ScalarAxis
bma = cifti.header.get_axis(1) # BrainModelAxis
for idx, (name, slc, bm) in enumerate(bma.iter_structures()): print((str(name), slc))
#('CIFTI_STRUCTURE_CORTEX_LEFT', slice(0, 29696, None))
#('CIFTI_STRUCTURE_CORTEX_RIGHT', slice(29696, None, None))
```
Label file:
```python
c = nib.load("Schaefer2018_100Parcels_7Networks_order.dlabel.nii")
c.shape	# (1, 59412)
c.header.get_axis(0) # LabelAxis
bma = c.header.get_axis(1) # BrainModelAxis
for idx, (name, slc, bm) in enumerate(bma.iter_structures()): print((str(name), slc))
# ('CIFTI_STRUCTURE_CORTEX_LEFT', slice(0, 32492, None))
# ('CIFTI_STRUCTURE_CORTEX_RIGHT', slice(32492, None, None))
c.header.get_axis(0).label #-> returns [{0 : ("???" , (1,1,1,0)), 1: ("area1",(R,G,B,A))}]
```



### Create a new Cifti (29k vertices describing left cortex excluding medial wall)

Each Cifti file needs an assocaited brain model (i.e. which areas are included in the cifti file, such as left hemisphere aka "LEFT_CORTEX", right hemisphere etc). This brain model we can take from an already existing file, such as the resting state run of an HCP subject:

```python
resting_state_run = "/data/hcp/sub-100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii"
img = nib.load(resting_state_run)

# extract brain model from file
# here we just take the left hemisphere. the brain models used in the HCP are special as they exclude the medial wall
# which reduces the number of greyordinates per hemisphere from 32492 to 29696

left_cortex_excluding_medial_wall = list(img.header.matrix._mims[1].brain_models)[0].vertex_indices._indices

# create a dummy array for the greyvoxel data, and only set the vertices that are not part of the medial wall to 1
mask = np.zeros((32492)); # the 32k resolution is standard; there is also higher ~ 164k?
np.put(mask, left_cortex_excluding_medial_wall, 1)
# create a new brain model with only the left cortex exluding medial wall
bm_left = nib.cifti2.BrainModelAxis.from_mask(mask, "LEFT_CORTEX")
```

Now we can finally create the actual image using this brain model. Multiple scalar images can be saved in a single cifti (i.e. curvature and cortical thickness, distrubution of different receptor types). Each of the scalar images needs a title:

```python
data = np.zeros((2, 29696)); # (n_scalars, n_vertices), for the current brainmodel n_vertices = 29696
# fill in / manipulate image data ...

dsnames = ["Curvature", "Thickness"];
cimg = nib.Cifti2Image(data, nib.cifti2.Cifti2Header.from_axes((nib.cifti2.ScalarAxis(dsnames), bm_left)))
# optionally save it to a file
cimg.to_filename("path/to/new/file.dscalar.nii(.gz)");
```

### Create a new scalar cifti (32k for left cortex, with medial wall, untested)
can be used on top of fsLR32k

```python
# create a dummy array for the greyvoxel data for the left hemisphere in 32k space, all set to one 1
mask = np.ones((32492));
# create a new brain model with only the left cortex
bm_left = nib.cifti2.BrainModelAxis.from_mask(mask, "LEFT_CORTEX")
data = np.zeros((2, 32492)); # (n_scalars, n_vertices), for the current brainmodel n_vertices = 32492
# fill in / manipulate image data ...
dsnames = ["Curvature", "Thickness"];
cimg = nib.Cifti2Image(data, nib.cifti2.Cifti2Header.from_axes((nib.cifti2.ScalarAxis(dsnames), bm_left)))
```

### Create a new scalar cifti (29k vertices for each left & right cortex)
medial wall is beeing excluded again, can be used on top of fsLR32k
```python
# mask out the medial wall on each side individually (using hcp-utils this time)
import hcp-utils as hcp
mask_l = np.zeros((32492)); np.put(mask_l, hcp.vertex_info['grayl'], 1)
mask_r = np.zeros((32492)); np.put(mask_r, hcp.vertex_info['grayr'], 1)

# and create the brain axes that exclude medial wall stuff
bma_left = nib.cifti2.BrainModelAxis.from_mask(mask_l, "LEFT_CORTEX")
bma_right = nib.cifti2.BrainModelAxis.from_mask(mask_r, "RIGHT_CORTEX")
bma_corticesLR = bma_left + bma_right
print(bma_corticesLR.nvertices) # {'CIFTI_STRUCTURE_CORTEX_LEFT': 32492, 'CIFTI_STRUCTURE_CORTEX_RIGHT': 32492}
for idx, (name, slc, bm) in enumerate(bma_corticesLR.iter_structures()):
    print((str(name), slc))
#('CIFTI_STRUCTURE_CORTEX_LEFT', slice(0, 29696, None))
#('CIFTI_STRUCTURE_CORTEX_RIGHT', slice(29696, None, None))
```

```python
caxes= (nib.cifti2.ScalarAxis(["MyelinBC_MSMAII"]), bma_corticesLR);
cheader=nib.cifti2.Cifti2Header.from_axes(caxes);
cimg = nib.Cifti2Image(myelin_data59k_LR , cheader)
```

### Create a new label cifti (29k vertices for each left & right cortex)

```python
# Create the Label axis with 3 'areas'
new_label_dict = {0 : ("unassigned" , (1,1,1,0)), 1: ("area1",(R,G,B,A)), 2: ("area2",(1,0.3,0.3,1)), ...}
lax = nib.cifti2.cifti2_axes.LabelAxis(["parcels"], new_label_dict)

# create the header and the cifti
cheader=nib.cifti2.Cifti2Header.from_axes((lax, bma_corticesLR));
# data_29kLR of shape (,59k) -> needs to be (1, 59k)
cimg = nib.Cifti2Image(np.expand_dims(data_29kLR, axis=0), cheader)
cimg.to_filename("data_29kLR.dlabel.nii");
```



### HCP: Load resting state run cifti-files and extract (left) cortex data

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

To get the resting state data from HCP files for the LH only, simply do:

```python
left_ctx_tseries29k; = fullbrain_tseries29k[:, 0: 0+29696];
```

Or more generally, for data in the HCP format (fsLR32k excluding medial wall, so actually having more like 29k), the `hcp-utils` toolbox can be used
```python
left_ctx_tseries = fullbrain_tseries[:, hcp.struct.cortex_left ];
# for non-HCP files that have 32k per hemisphere, do:
run29k_L = run32k_L[:, hcp.vertex_info['grayl'] ]
# if both left and right are concatenated in the run (2x32k = 64k), one can do the following:
cortexLR = list(hcp.vertex_info['grayl']) + list(hcp.vertex_info['grayr'] + 32492)
run59k_LR = run64k_L[:, hcp.vertex_info['grayl'] ]
```

For displaying on 32k mesh files of i.e. the cortical surface, we need to project onto the relevant vertices in this 32k space (others are left at 0):
```python
# to project back onto 32k 
left_cortex_excluding_medial_wall = list(img.header.matrix._mims[1].brain_models)[0].vertex_indices._indices
left_ctx_tseries32k = np.zeros((32492));
left_ctx_tseries32k[left_cortex_excluding_medial_wall] = left_ctx_tseries29k;

# to project back onto 32k using hcp-utils
left_ctx_tseries32k = hcp.left_cortex_data(left_ctx_tseries29k;)    # returns 32k version
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


### Get outline of a surface ROI

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

### Combine Gifti into Cifti

```shell

#works similiarly both labels and func: (R.label.gii, L.label.gii) -> (dlabel.nii) and (R.func.gii, L.func.gii) -> (dscalar.nii)
# wb_command -cifti-create-dense-scalar <output> [-part <part.gifti>]
wb_command -cifti-create-dense-scalar juelich_atlas_v29.maxprob.enc.32kfslr.LR.dscalar.nii -left-metric juelich_atlas_v29.maxprob.enc.32kfslr.L.func.gii -right-metric juelich_atlas_v29.enc.32kfslr.R.func.gii

# wb_command -cifti-create-label <output> [-part <part.gifti>]
wb_command -cifti-create-label juelich_atlas_v29.maxprob.enc.32kfslr.LR.dlabel.nii -left-label juelich_atlas_v29.maxprob.enc.32kfslr.L.label.gii -right-label juelich_atlas_v29.enc.32kfslr.R.label.gii

# ... and timeseries: -cifti-create-dense-timeseries
```
* create label needs a label mapping (`gimg.labeltable.labels`) in both of the metric (`*.label.gii`) files, as only those parcels are included in the combined image that are described in the mapping

## Gifti


* `.surf.gii` files (also called `Geometry` or `surface` files) provide the 3D brain mesh/object backbone. Contains the mapping from greyordinate vertex to its location in 3D. They are also referred to as 
* `shape.gii` or `func.gii` (also called `metric` files) contain actual data, such as cortical thickness or functional activations etc. (similiar to what can be contained in cifti files?). func is usually used for multi-timepoint data 
* `label.gii` contain integer data, indicating the parcel each vertex belongs to; + a mapping from these integer values to parcel labels (a name and a color)


Notably, other software may put data arrays (the equivalent of a metric file) into the same file as the geometry information. The connectome workbench doesnt support this.

For more info, check the posts by [Emma Robinson on Gifti/HCP surface files](https://emmarobinson01.com/2016/02/10/unofficial-guide-to-the-hcp-surface-file-formats/) and both posts by Joset (Jo) A. Etzel on [NIfTI, CIFTI, GIFTI in the HCP and Workbench](http://mvpa.blogspot.com/2014/03/nifti-cifti-gifti-in-hcp-and-workbench.html) and on [Conversion from Volumetric to Surface](https://mvpa.blogspot.com/2018/02/connectome-workbench-making-surface.html).




### Convert annoations between surface data types

*From gifti to cifti (within FS_LR32k)* 
```shell
# use cifti-create functions: -cifti-create-dense-timeseries, -cifti-create-label
# see also last sections of CIFTI (Separate Cifti files into gifti / Combine Gifti into Cifti)
wb_command -cifti-create-dense-scalar juelich_atlas_v29.maxprob.enc.32kfslr.LR.dscalar.nii -left-label juelich_atlas_v29.maxprob.enc.32kfslr.L.func.gii -right-label juelich_atlas_v29.enc.32kfslr.R.func.gii
```

*From cifti to gifti (within FS_LR32k)* 
```shell
wb_command -cifti-separate Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_LEFT Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii
```

*From FS_LR32k gifti to fsaverage164 gifti* [[fs_L-to-fs_LR template]](https://github.com/Washington-University/HCPpipelines/blob/master/global/templates/standard_mesh_atlases/fs_L/fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii "[fs_L-to-fs_LR template]")
```shell
wb_command -label-resample <label-in> <current-sphere> <new-sphere> <method> <label-out>
wb_command -label-resample Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii L.sphere.32k_fs_LR.surf.gii fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii BARYCENTRIC left.fsaverage164.label.gii

# "The ADAP_BARY_AREA method is recommended for label data, because it should be better at resolving vertices that are near multiple labels, or in case of downsampling"
wb_command -label-resample Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii L.sphere.32k_fs_LR.surf.gii fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii ADAP_BARY_AREA left.fsaverage164.label.gii
```

Related functions are [`-surface-resample`](https://www.humanconnectome.org/software/workbench-command/-surface-resample) and [`-metric-resample`](https://www.humanconnectome.org/software/workbench-command/-metric-resample) for the two other respective gifti file types.

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


### Add annotations to a func.gii (transform func.gii -> label.gii)

```python
# load unannotated (functional) gifti (i.e. resultung from a conversion from volmMNI -> fsLR32k gifti)
gimg = nib.load(r"JULICH_BRAIN_CYTO_29_MNI152_2009C_nA.maxprob.enc.32kfslr.L.func.gii")

# look at the already existing labels:
len(gimg.labeltable.labels) #->1 (labels is a list containing only one element here)
gimg.labeltable.labels[0].to_xml()
#> <Label Key="0" Red="1.0" Green="1.0" Blue="1.0" Alpha="0.0">???</Label>
# this basically indicates all voxels with value 0 are not assigned to any area

# labels for indices are extracted from somewhere or manually assigned:
# indices,labels = (list<int>, list<string>) 

for (idx,label) in zip(indices,labels): print(str(idx) +": " + label)
"""
1: Frontal-to-Temporal-I (GapMap) left
2: Frontal-to-Temporal-I (GapMap) right
3: Ch 123 (Basal Forebrain) left
4: Ch 123 (Basal Forebrain) right
"""

for (idx,label) in zip(indices,labels):
  r,g,b = list(np.random.rand(3).round(2))		 # create some randome color for the annotation
  l = nib.gifti.gifti.GiftiLabel(idx, r,g,b,1)			# create a GiftiLabel for the index (i.e. 3)
  l.label = label;												 # actually assign a name/label to the parcel, i.e. "Ch 123 (Basal Forebrain) left"
  gimg.labeltable.labels = gimg.labeltable.labels + [l]		# extend the existing labels in for the image

# now save the gifti as .label.gii
nib.save(gimg,"JULICH_BRAIN_CYTO_29_MNI152_2009C_nA.maxprob.enc.32kfslr.L.label.gii") 
```


### Create a gifti (.label.gii) from scratch (or using an annot file)

```python
import nibabel as nib

annotfn   = r"\lausanne2018_fsaverage\lh.lausanne2008.scale1.annot"  
giifn     = r"\lausanne2018_fsaverage\lh.lausanne2008.scale1.gii"  
label_out = r"lausanne08.sc1.fslr32k.L.label.gii"  

## load annot file and convert directly to gifti  
fslabels, fsctab, fsnames = nib.freesurfer.io.read_annot(annotfn)  
#fslabels: 163842 labeled voxels in L, i.e.: [10 18  8  7 10 17 23 20  0  1] (35 unique labels)
#fsctab if of length 35, with fsctab[0] beeing i.e. [250 250 250 0 16448250] (RGBT+)
#fsnames (len of 35), i.e.: [b'unknown', b'lateralorbitofrontal', b'parsorbitalis']

# usually a gifti describes a single anatomical strcture
# this has to be given as metadata, as otherwise i.e. the HCP workbench doesnt know how to place the data onto a mesh 
gmeta =  nib.gifti.GiftiMetaData(nib.gifti.GiftiNVPairs("AnatomicalStructurePrimary", "CortexLeft"))  
# data should be the same size and order of the vertices contained within the relevant surface mesh (here: the annot file described the left hemisphere of fasaverage(7) which has 164k vertices (163842)
ga = nib.gifti.GiftiDataArray(data=fslabels, intent="NIFTI_INTENT_NORMAL", datatype="NIFTI_TYPE_INT32");  # dtpe can also be *_FLOAT32 (=default?)

# for .label.gii files we should also provide a lable table (this can also aparrently contain unassigned labels; these are merged when a cifti is created from two hemisphere giftis
lt = nib.gifti.GiftiLabelTable();  
for i in range(len(fsnames)):  
  r,g,b = np.array(fsctab[i, :3])/255	# colors should be in range [0, 1]  
  l = nib.gifti.gifti.GiftiLabel(i, r,g,b,1) # create a GiftiLabel for the index (i.e. 3)  
  l.label = fsnames[i].decode("utf-8"); # actually assign a name/label to the parcel, i.e. "Ch 123 (Basal Forebrain) left"  
  lt.labels = lt.labels + [l]    # extend the existing labels in for the image  
  
# finally create and save the gifti
g = nib.gifti.GiftiImage(labeltable=lt, darrays=[ga], meta=gmeta)  
g.to_filename(giifn)
```

```python
# do the other hemisphere too ...

# optionally join the two resulting fsLR32k .gii files again to a joint cifti (.dlabel.nii)  
labelL = os.path.join(adir, "fsLR32k", f"lausanne08.sc{scale}.fslr32k.L.label.gii")  
labelR = os.path.join(adir, "fsLR32k", f"lausanne08.sc{scale}.fslr32k.R.label.gii")  
cmd = f"wb_command -cifti-create-label {out_cifti} -left-label {labelL} -right-label {labelR}"

```

### Get an ROI from a Gifti parcellation

In the HCP data a parcellation is i.e. given for each subject under:
`hcp-subjectname/MNINonLinear/fsaverage_LR32k/subjectname.L.aparc.32k_fs_LR.label.gii`
(here the one in the FS_LR32k space is used)

```python
import nibabel as nib

AnatLabels = nib.load(r"/data/hcp/hcp-subjectname/MNINonLinear/fsaverage_LR32k/subjectname.L.aparc.32k_fs_LR.label.gii")
AnatLabelsData= AnatLabels.darrays[0].data
# -> array([10, 29, 24, ..., 15, 15, 15] of shape (32492,)

# the mapping of those values to labels is given in 
AnatLabels.get_labeltable().labels

# if we were interested in the IFG pars triangularis, we would use the vertices with value 20:
AnatLabels.get_labeltable().labels[20].label 
#-> u'L_parstriangularis'

triangularis_mask = AnatLabelsData == 20;
```




## Misc: working with surface data in a 2D plane


```python


```

<br/><br/>

# Toolboxes

## hcp-utils

https://rmldj.github.io/hcp-utils/

```python
!pip install hcp_utils
import hcp_utils as hcp
img = nib.load('path/to/fMRI_data_file.dtseries.nii')
X = img.get_fdata()
X.shape     # e.g. (700, 91282)
X_hipL = X[:, hcp.struct.hippocampus_left]
X_hipL.shape    # (700, 764)

Xp = hcp.parcellate(Xn, hcp.yeo7)
Xp.shape    # (700, 7)
plotting.view_surf(mesh_sub.inflated, 
    hcp.cortex_data(hcp.unparcellate(Xp[29], hcp.yeo7)), 
    threshold=0.1, bg_map=mesh_sub.sulc)
```

```python
import hcp_utils as hcp
# the HCP file s1200_sulc i.e. contains 29k for left cortex and 29k for right cortex
s1200_sulc_L = s1200_sulc[0,hcp.struct.cortex_left]     # returns 29k version (LH excluding medial wall)
s1200_sulc_L = hcp.left_cortex_data(s1200_sulc[0,:])    # returns 32k version
#s1200_sulc_L = hcp.left_cortex_data(s1200_sulc[0,hcp.struct.cortex_left]) # returns 32k version
```

```python
# to project back onto 32k using hcp-utils
left_ctx_tseries32k = hcp.left_cortex_data(left_ctx_tseries29k;)    # returns 32k version
# for non-HCP files that have 32k per hemisphere, do:
run29k_L = run32k_L[:, hcp.vertex_info['grayl'] ]
# if both left and right are concatenated in the run (2x32k = 64k), one can do the following:
cortexLR = list(hcp.vertex_info['grayl']) + list(hcp.vertex_info['grayr'] + 32492)
run59k_LR = run64k_L[:, hcp.vertex_info['grayl'] ]
```


## brainspace

https://github.com/MICA-MNI/BrainSpace

https://brainspace.readthedocs.io/en/latest/

```python
syt20 = nib.load("Schaefer2018_400Parcels_7Networks_order_Tian_Subcortex_S2.dlabel.nii")
syt20_LR29k = syt20.get_fdata()[0,hcp.struct.cortex] # shape 59412 ~ 2*29k

from brainspace.utils.parcellation import reduce_by_labels, map_to_labels
print("Number of parcels in ctx atlas: ", len(np.unique(syt20_LR29k)))
#Number of parcels in ctx atlas:  401
mask = syt20_LR29k!=0; # only where the atlas is not zero
timeseries_red = reduce_by_labels(timeseries[mask], syt20_LR29k[mask], axis=1, red_op='mean') #resutling shape: (400, 4800)
grad=map_to_labels(timeseries_red, syt20_LR29k, mask=mask, fill=np.nan)
```

## netneurotools

This toolbox is a collection of functions written in Python that get frequent usage in the Network Neuroscience Lab, housed in the Brain Imaging Centre at McGill University.

https://github.com/netneurolab/netneurotools/

https://netneurotools.readthedocs.io/en/latest/generated/netneurotools.datasets.fetch_cammoun2012.html

https://github.com/netneurolab/netneurotools/blob/master/resources/generate_atl-cammoun2012_surface.py

```python
from netneurotools.datasets import *
fetch_cammoun2012("fslr32k", data_dir=r"C:\tmp\OwnCloud\res\atlases\lausanne", verbose=0)
```

## neuromaps

The neuromaps toolbox is designed to help researchers make easy, statistically-rigorous comparisons between brain maps (or brain annotations).

https://github.com/netneurolab/neuromaps


```python
from neuromaps.datasets import available_annotations
for annotation in available_annotations(): print(annotation)
# ('laurikainen2018', 'fmpepd2', 'MNI152', '1mm') <source>, <description>, <space>, <resolution>
annotation = fetch_annotation(source='neurosynth') # tags=
```

## AbaGen

https://github.com/rmarkello/abagen

In 2013, the Allen Institute for Brain Science released the Allen Human Brain Atlas, a dataset containing microarray expression data collected from six human brains (Hawrylycz et al., 2012). [...] The current Python package, abagen, aims to provide reproducible workflows for processing and preparing the AHBA microarray expression data for analysis.



<br/><br/><br/><br/><br/><br/><br/><br/><br/>

# FSaverage spaces

Information taken from [this mailing list thread](https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/msg28867.html):

> Hi Ludovico, I'm not sure exactly what your pipeline is, but you should use  fsaverage for vol2surf, then go from fsaverage to fsaverage4 or 5 with  mri_surf2surf. The reason is that if you go directly to 4 or 5, you will miss data between the vertices. The distance between the vertices for **fsaverage4 is about 5.5mm** and it is about **3mm for fsaverage5**, **1.4 for fsaverage6**, and **.7mm for fsaverage** (these are all for the white surface).

So it seems as if the standard fsaverge could be considered fsaverage7.

>  Which surface is "6mm" defined on? Is it .sphere, .white, or something else?

> It is based on the white. Notice that since you are using an average subject (fsaverage6) there is some scaling that happens because an average subject will have a surface area that is much less than the average of the subjects that went into. The average of the subjects is kept when the average subject is created and used when the number of iterations is computed.

**fsaverage4** 
seems to have 2329 vertices in one hemisphere.

**fsaverage5**
Data are defined on a surface with 10,024 vertices (FreeSurfer **fsaverage 5**). One shows the standard averaging referred to as Mean, the averaging after Gaussian smoothing is referred to as Mean (S) (mean after Gaussian smoothing with

**fsaverage164**

probably corresponds to fsaverage == fsaverage7; model calculation:

* fsaverge5 if it was squared would have $\sqrt(10024)$ ~ 100 vertices on each side, with between vertex distance of 3mm. Hence each side would roughly of length 100vertices x 3mm = 300mm
* we know for fsaverage(7) the mean vertex distance is 0.7mm, that is each side of length 300mm would contain ~ 428 vertices (300mm/0.7mm/vertex). So the total number of vertices can be roughly estimated as $428^2 = 183184$ which is close to 164k


**info from nilearn**
‘fsaverage3’: the low-resolution fsaverage3 mesh (642 nodes)
‘fsaverage4’: the low-resolution fsaverage4 mesh (2562 nodes)
‘fsaverage5’: the low-resolution fsaverage5 mesh (10242 nodes)
‘fsaverage6’: the medium-resolution fsaverage6 mesh (40962 nodes)
‘fsaverage7’: same as ‘fsaverage’
‘fsaverage’: the high-resolution fsaverage mesh (163842 nodes)
https://nilearn.github.io/modules/generated/nilearn.datasets.fetch_surf_fsaverage.html#nilearn.datasets.fetch_surf_fsaverage




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


## HCP MMP1 on fsaverge

https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446

## HCP MMP1 in AFNI

FNI will now include the Human Connectome Project brain atlas in MNI space, originally from Glasser, et al, Nature 2016. It will be available by default in AFNI's whereami output, both in the GUI and on the command line. 

https://openwetware.org/wiki/Beauchamp:CorticalSurfaceHCP
https://afni.nimh.nih.gov/pub/dist/atlases/MNI_HCP/

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



