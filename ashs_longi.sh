#!/bin/bash
set -x -e

# Inputs: directories where baseline and followup ASHS have been run
ID=${1?}
ASHSBL=${2?}
ASHSFU=${3?}
SIDE=${4?}
ATLAS=${5?}
WORK=${6?}

# Baseline and followup scans
MPRAGE_BL=$ASHSBL/mprage.nii.gz
MPRAGE_FU=$ASHSFU/mprage.nii.gz

# TSE - fine to use chunks
TSE_BL_RAW=$ASHSBL/tse_native_chunk_${SIDE}.nii.gz
TSE_FU_RAW=$ASHSFU/tse_native_chunk_${SIDE}.nii.gz

# T1-T2 matrices
MAT_T1T2_BL=$ASHSBL/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat
MAT_T1T2_FU=$ASHSFU/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat

# ASHS segmentation
ASHS_SEG_BL=$ASHSBL/bootstrap/fusion/lfseg_heur_${SIDE}.nii.gz
ASHS_SEG_FU=$ASHSFU/bootstrap/fusion/lfseg_heur_${SIDE}.nii.gz

# Intermediates
MAT_MPRAGE_BLFU=$WORK/mprage_whole_rigid.mat
MAT_TSE_BLFU_INIT=$WORK/tse_rigid_init.mat
MAT_TSE_BLFU_NATIVE=$WORK/tse_rigid_native.mat
MASK_BL_NATIVE=$WORK/bl_mask.nii.gz
MASK_BL_HW=$WORK/bl_mask_hw.nii.gz
MASK_BL_HW_CROP=$WORK/bl_mask_hw_crop.nii.gz
MASK_BL_HW_CROP_IT2=$WORK/bl_mask_hw_crop_it2.nii.gz
REFSPACE=$WORK/refspace.nii.gz

MAT_M1=$WORK/M1.mat
MAT_M2=$WORK/M2.mat
MAT_B=$WORK/B.mat
MAT_B_SQRT=$WORK/B_sqrt.mat
MAT_B_SQRT_INV=$WORK/B_sqrt_inv.mat
MAT_Z1=$WORK/Z1.mat
MAT_Z2=$WORK/Z2.mat
MAT_Y1=$WORK/Y1.mat
MAT_Y2=$WORK/Y2.mat

MAT_Z2_Y1_INV=$WORK/Z2_Y1_inv.mat
MAT_Y2_Z1_INV=$WORK/Y1_Z2_inv.mat
MAT_B_SQRT_IT2=$WORK/B_sqrt_it2.mat
MAT_B_SQRT_INV_IT2=$WORK/B_sqrt_inv_it2.mat
MAT_Z1_IT2=$WORK/Z1_it2.mat
MAT_Z2_IT2=$WORK/Z2_it2.mat

TSE_BL=$WORK/tse_bl.nii.gz
TSE_FU=$WORK/tse_fu.nii.gz
TSE_BL_HW=$WORK/reslice_tse_bl_to_hw.nii.gz
TSE_FU_HW=$WORK/reslice_tse_fu_to_hw.nii.gz
TSE_BL_HW_IT2=$WORK/reslice_tse_bl_to_hw_it2.nii.gz
TSE_FU_HW_IT2=$WORK/reslice_tse_fu_to_hw_it2.nii.gz

MASK_JLF=$WORK/mask_jlf.nii.gz

# Make work directory
mkdir -p $WORK

<<'SKIP'


SKIP

# Resample the TSE images and organize as RSA
c3d $TSE_BL_RAW -resample 100x100x500% -swapdim RSA -o $TSE_BL
c3d $TSE_FU_RAW -resample 100x100x500% -swapdim RSA -o $TSE_FU

# Perform whole-brain T1 registration between the T1 scans without a mask
greedy -d 3 -a -dof 6 \
  -i $MPRAGE_BL $MPRAGE_FU -m NCC 4x4x4 -n 40x80x40x0 \
  -o $MAT_MPRAGE_BLFU -ia-image-centers

# Create a registration mask from the ASHS output (baseline is ok)
c3d $TSE_BL $ASHS_SEG_BL -thresh 1 inf 1 0 -dilate 1 8x8x0 \
  -reslice-identity -o $MASK_BL_NATIVE

# Peform mask-based registration between the T2 scans using the T2->T1 transforms
# from ASHS as the starting point
#
# TODO: should the mask be cropped to account for partial overlap?
c3d_affine_tool $MAT_T1T2_BL -inv $MAT_MPRAGE_BLFU -mult $MAT_T1T2_FU -mult -o $MAT_TSE_BLFU_INIT
greedy -d 3 -a -dof 6 -m NCC 3x3x1 \
  -i $TSE_BL $TSE_FU -gm $MASK_BL_NATIVE -o $MAT_TSE_BLFU_NATIVE -ia $MAT_TSE_BLFU_INIT

# Create reference space
c3d $TSE_BL -origin 0x0x0mm -orient RSA -o $REFSPACE

# Matrices from reference space to T2 space
c3d_affine_tool -sform $TSE_BL -sform $REFSPACE -inv -mult -info -o $MAT_M1
c3d_affine_tool -sform $TSE_FU -sform $REFSPACE -inv -mult -info -o $MAT_M2

# B matrix
c3d_affine_tool $MAT_M2 -inv $MAT_TSE_BLFU_NATIVE -mult $MAT_M1 -mult -o $MAT_B \
  -sqrt -o $MAT_B_SQRT -inv -o $MAT_B_SQRT_INV

# Matrices that apply half-way transformations
c3d_affine_tool $MAT_M1 $MAT_B_SQRT_INV -mult -o $MAT_Z1
c3d_affine_tool $MAT_M2 $MAT_B_SQRT -mult -o $MAT_Z2


# Test that the images line up
mkdir -p $WORK/test
c3d $REFSPACE $TSE_BL -reslice-matrix $MAT_Z1 -o $TSE_BL_HW
c3d $REFSPACE $TSE_FU -reslice-matrix $MAT_Z2 -o $TSE_FU_HW
c3d $REFSPACE $MASK_BL_NATIVE -int 0 -reslice-matrix $MAT_Z1 -o $MASK_BL_HW

# Compute a cropped mask
c3d $REFSPACE -as R $TSE_BL -thresh -inf inf 1 0 \
  -pad 0x0x1vox 0x0x1vox 0 -dilate 0 1x1x1 -reslice-matrix $MAT_Z1 -popas K1 \
  -push R $TSE_FU  -thresh -inf inf 1 0 \
  -pad 0x0x1vox 0x0x1vox 0 -dilate 0 1x1x1 -reslice-matrix $MAT_Z2 -popas K2 \
  -push K1 -push K2 -times $MASK_BL_HW -times -o $MASK_BL_HW_CROP

# Perform a registration between halfway images - to account for residual misreg
greedy -d 3 -dof 6 \
  -i $TSE_BL_HW $TSE_FU -gm $MASK_BL_HW_CROP \
  -o $MAT_Y2 -ia $MAT_Z2 -a -m NCC 4x4x4 -n 40x40

greedy -d 3 -dof 6 \
  -i $TSE_FU_HW $TSE_BL -gm $MASK_BL_HW_CROP \
  -o $MAT_Y1 -ia $MAT_Z1 -a -m NCC 4x4x4 -n 40x40

# Compute the new Z (native to halfway reference) matrices
c3d_affine_tool $MAT_Z2 $MAT_Y1 -inv -mult -info -o $MAT_Z2_Y1_INV
c3d_affine_tool $MAT_Y1 $MAT_Z1 -inv -mult -info -o $MAT_Y2_Z1_INV
c3d_affine_tool $MAT_M2 -inv $MAT_Z2_Y1_INV -mult $MAT_M1 -mult \
  -sqrt -o $MAT_B_SQRT_IT2 -inv -o $MAT_B_SQRT_INV_IT2 
c3d_affine_tool $MAT_M1 $MAT_B_SQRT_INV_IT2 -mult -o $MAT_Z1_IT2
c3d_affine_tool $MAT_M2 $MAT_B_SQRT_IT2 -mult -o $MAT_Z2_IT2

# Reslice to reference space again
c3d $REFSPACE $TSE_BL -reslice-matrix $MAT_Z1_IT2 -o $TSE_BL_HW_IT2
c3d $REFSPACE $TSE_FU -reslice-matrix $MAT_Z2_IT2 -o $TSE_FU_HW_IT2

# Compute a cropped mask
c3d $REFSPACE -as R $TSE_BL -thresh -inf inf 1 0 \
  -pad 0x0x1vox 0x0x1vox 0 -dilate 0 1x1x1 -reslice-matrix $MAT_Z1 -popas K1 \
  -push R $TSE_FU  -thresh -inf inf 1 0 \
  -pad 0x0x1vox 0x0x1vox 0 -dilate 0 1x1x1 -reslice-matrix $MAT_Z2 -popas K2 \
  -push K1 -push K2 -times \
  -push R $MASK_BL_NATIVE -int 0 -reslice-matrix $MAT_Z1_IT2 \
  -times -o $MASK_BL_HW_CROP_IT2

# Perform ASHS voting
mkdir -p $WORK/watlas
for tid in $(ls $ATLAS/train); do

  # Bootstrap directory for BL and FU
  BOOT_BL=$ASHSBL/bootstrap/tseg_${SIDE}_${tid}
  BOOT_FU=$ASHSFU/bootstrap/tseg_${SIDE}_${tid}

  # Get the chain of transforms that takes the atlas into the subject T2 space
  CHAIN_BL="$MAT_Z1_IT2 $BOOT_BL/sqrt_inv.mat,-1 $BOOT_BL/greedy_warp.nii.gz $BOOT_BL/sqrt_fwd.mat"
  CHAIN_FU="$MAT_Z2_IT2 $BOOT_FU/sqrt_inv.mat,-1 $BOOT_FU/greedy_warp.nii.gz $BOOT_FU/sqrt_fwd.mat"

  # Apply the chains to the MRI and labels
  greedy -d 3 \
    -rf $REFSPACE \
    -rm $ATLAS/train/${tid}/tse.nii.gz $WORK/watlas/atlas_bl_${tid}.nii.gz \
    -ri LABEL 0.2vox \
    -rm $ATLAS/train/${tid}/seg_${SIDE}.nii.gz $WORK/watlas/atseg_bl_${tid}.nii.gz \
    -r $CHAIN_BL

  greedy -d 3 \
    -rf $REFSPACE \
    -rm $ATLAS/train/${tid}/tse.nii.gz $WORK/watlas/atlas_fu_${tid}.nii.gz \
    -ri LABEL 0.2vox \
    -rm $ATLAS/train/${tid}/seg_${SIDE}.nii.gz $WORK/watlas/atseg_fu_${tid}.nii.gz \
    -r $CHAIN_FU

done


# Compute the mask for LF more carefully, excluding edge regions of both segmentations
# where there is more uncertainty
c3d \
  $REFSPACE -as R $ASHS_SEG_BL -thresh 1 inf 1 0 -as X -dup \
  -pad 0x0x1 0x0x1 0 -dilate 1 300x300x0 -dilate 0 0x0x1 \
  -int 0 -reslice-identity -push X \
  -dilate 1 4x4x0 -times \
  -reslice-matrix $MAT_Z1_IT2 \
  -push R $ASHS_SEG_FU -thresh 1 inf 1 0 -as X -dup \
  -pad 0x0x1 0x0x1 0 -dilate 1 300x300x0 -dilate 0 0x0x1 \
  -int 0 -reslice-identity -push X \
  -dilate 1 4x4x0 -times \
  -reslice-matrix $MAT_Z2_IT2 \
  -times -o $MASK_JLF

# Directory for the posteriors
mkdir -p $WORK/lfpost

# Do the label fusion
~/tk/malf/xc64rel/label_fusion 3 \
  -g $WORK/watlas/atlas_*.nii.gz \
  -l $WORK/watlas/atseg_*.nii.gz \
  -rp 3x3x3 -rs 3x3x3 -M $MASK_JLF \
  -p $WORK/lfpost/post_bl_%03d.nii.gz \
  $TSE_BL_HW_IT2 $WORK/ashsseg_bl_hw.nii.gz

~/tk/malf/xc64rel/label_fusion 3 \
  -g $WORK/watlas/atlas_*.nii.gz \
  -l $WORK/watlas/atseg_*.nii.gz \
  -rp 3x3x3 -rs 3x3x3 -M $MASK_JLF \
  -p $WORK/lfpost/post_fu_%03d.nii.gz \
  $TSE_FU_HW_IT2 $WORK/ashsseg_fu_hw.nii.gz


# Use MRF to create regularized segmentations
mkdir -p $WORK/mrf
for i in 000 001 002 003 004 006 007 008 009 011 012 013; do 
  
  # Make sure all the posteriors exist
  for what in bl fu; do
    if [[ ! -f $WORK/lfpost/post_${what}_${i}.nii.gz ]]; then
      c3d $REFSPACE -scale 0 -o $WORK/lfpost/post_${what}_${i}.nii.gz 
    fi
  done

  # Create the 4D posteriors
  c4d $WORK/lfpost/post_bl_$i.nii.gz $WORK/lfpost/post_fu_$i.nii.gz \
    -tile w -o $WORK/mrf/post_blfu_${i}.nii.gz

done

# Create a 4D mask
c4d $MASK_JLF -dup -tile w -o $WORK/mrf/mask_mrf.nii.gz

# Perform MRF with different lambdas
for lam in 0.01 0.02 0.04 0.06 0.08 1.0; do

  lamcode=$(echo $lam | awk '{printf "lam_%04d",1000 * $1}')

  c4d \
    $WORK/mrf/post_blfu_*.nii.gz $WORK/mrf/mask_mrf.nii.gz -popas M \
    -foreach -clip 0 1 -push M -replace 0 NaN -times -endfor \
    -vote-mrf VA $lam -o $WORK/mrf/vote4d_${lamcode}.nii.gz \
    -replace 1 0 2 1 3 2 4 3 5 4 6 5 7 7 8 8 9 9 10 11 11 12 12 13 \
    -slice w 0:-1 -oo $WORK/mrf/vote_bl_${lamcode}.nii.gz $WORK/mrf/vote_fu_${lamcode}.nii.gz

done






