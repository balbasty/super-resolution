# Denoising / Super-Resolution

This code denoises or super-resolves multi-contrast MR images.

It inverts a generative model of the data with a multi-channel 
total-variation prior. It handles multiple contrasts with different 
resolutions and field-of-views. It can also take advantages of multiple 
repeats of a same contrast.

## Dependencies

This code is written in Matlab and depends on the **development branch**
of the [SPM12](https://www.fil.ion.ucl.ac.uk/spm/) software, which 
should be on your Matlab path. Access to this branch is usually 
restricted to people who work at the 
[FIL](https://www.fil.ion.ucl.ac.uk/spm/local/). If you don't but would 
like to try this software, please send us an email to get a copy.

For increased speed, it is advised to recompile SPM with OpenMP 
activated. This can be done by following 
[these instructions](https://en.wikibooks.org/wiki/SPM), except that the 
option `USE_OPENMP=1` must be specified. For example, on linux:
```{shell}
cd /home/login/spm12/src
make distclean
make USE_OPENMP=1 && make install
make external-distclean
make external && make external-install
```

## Usage

### SPM Batch

If the super-resolution is copied into the `toolbox/` folder of SPM, it
can be used within the SPM batch system, which can be accessed by clicking 
the **Batch** Button in the main window. The super-resolution tool can 
then be loaded from  `SPM > Tools > Denoising / Super-Resolution`.

The top-level field allows the user to chooose between the super-resolution 
and denoising modes. Most other options are identical between these two 
modes, but the super-resolution lets the user specify a few more 
parameters, such as the target voxel-size and the slice profiles.

### Graphical

The tool can also be run using a small graphical interface outside of SPM.
Assuming that you have a series of noisy images in nifti format, 
and that you like graphical interfaces, the most simple usage is:
```{matlab}
>> sr_fit;
```
You will be asked:
- **Number of contrasts**: the number of different MR contrasts (T1w, T2w, ...)
- **Mode** (*denoising*|*super-resolution*): super-resolution allows 
  MRIs to be upsampled to a finer grid (*i.e.*, smaller voxel size).
  Note that super-resolution *implies* denoising.
  - **Voxel Size**: Target voxel size (only for *super-resolution*).
  - **Thick slice** (*yes*|*no*): The thick-slice mode can be used for 
    anisotropic acquisitions, with high in-plane and low through-plane 
    resolution (only for *super-resolution*).
- **Regularisation** (*5E3*): Higher = stronger denoising
- **Number of iterations** (*10*): Higher = closer to the optimum

You will then be asked to select all input files for each contrast. 
You can select multiple files (*i.e.*, repeats) for each contrast.
Finally, you will be asked to choose an output directory.

The algorithm will then run and write the results on disk. 
There will be as many denoised/super-resolved images as contrasts.

### Command line

If more flexibility is required, `sr_fit` can be used as a function:
```{matlab}
>> [out,in] = sr_fit(dat, opt);
```
The input `dat` can be either:
- A cell of cells of nifti filenames. The outer loop corresponds to 
  contrasts and the inner loop corresponds to repeats of a given 
  contrast.
- A 4D numeric or file array. This works in the denoising case only, 
  when all inputs are registered and reslices to the same lattice.
- Empty. In this case, files can be selected through the GUI.

The outputs are:
- `out` is a structure with fields:
  - `dat` - a file array containing the denoised images
  - `rls` - a file array containing the L1 weights
  - `mat` - the output orientiation matrix
  - `dim` - the output dimensions
  - `lam` - the corrected regularisation factor for each contrast).
- `in` is a cell of cells of structures with fields:
  - `dat` - a file array containing the observed data
  - `mat` - the input orientation matrix
  - `dim` - the input dimensions
  - `var` - the observation uncertainty (based on the data type)
  - `lam` - the estimated noise precision
  - `mu`  - the estimated mean tissue intensity 
(used to corrected the regularisation factors)
- If no output argument is asked (`>> r_fit(dat, opt);`), an output 
  folder is asked for through the GUI and the output images are writen 
  on disk.

Several options are available, as fields of the `opt` structure:
| Name         | Range | Default | Description |
|--------------|-------|---------|-------------|
| `mode`       | 'denoise' / 'superres' | 'denoise' | |
| `log`        | bool       | false | Solve for-log images (makes the prior scale-indepdendent) |
| `itermax`    | int > 0    | 10 | Maximum number of iterations |
| `tolerance`  | float > 0  | 1E-3 | Gain threshold for early stopping |
| `reg.mode`   | 0 / 1 / 2  | 1 | Regularisation mode: 0=None/1=L1/2=L2 |
| `reg.value`  | float > 0  | 5E4 | Regularisation factor per channel |
| `reg.smo`    | float > 0  | 1E-3 | RLS smoothing term |
| `vs`         | float > 0  | NaN | Target voxel size. By default, it is the average input voxel size. |
| `fov`        | int >= 0   | 0 | Target field of view. If 0, an average orientation matrix is computed from the input matrices. If N>0, the orientation matrix of the N-th input image is used. |
| `coreg.do`   | bool       | true | Start by coregistering all input volumes. |
| `coreg.fwhm` | float > 0  | [21 14 7] | Series of FWHM used to smooth the joint histogram |
| `threads`    | int >= 0 / 'matlab' / 'automatic' |  'matlab' | Number of threads used by both Matlab and SPM: 'matlab' = Matlab's current settings / 'automatic = Matlab's automatic setting |
| `verbose`    | int >= 0   | 1 | Verbosity level: 0=quiet / 1=print / 2=plot |
| `out.mem`    | 'map' / 'load' | 'map' | Memory map output data (slower but saves RAM) |
| `out.folder` | char       | '.' | Output folder. |
| `slice.dir`  | 'thickest' / 'all' | 'thickest' | Which directions are 'slice-selection' directions? |
| `slice.gap`  | float >= 0 | 1/3 | Gap between slices in the slice direction(s) |
| `armijo`     | float >= 0 | [2 1] | Series of damping factors used for Gauss-Newton . |
| `solver.fmg`   | int >= 0 | [2 2] | Number of cycles (1) and relaxation iterations (2) used by the multigrid solver . |
| `solver.cg`    | int >= 0 | 50    | Number of conjugate gradient iterations . |
| `solver.relax` | int >= 0 | 0     | Series of relaxation iterations . |

## References

The use of Multi-Channel Total-Variation as a prior for MR 
super-resolution is described in:

- **MRI Super-Resolution using Multi-Channel Total Variation.**  
[Mikael Brudfors](mailto:brudfors@gmail.com), [Yaël Balbastre](mailto:y.balbastre@ucl.ac.uk), [Parashkev Nachev](mailto:p.nachev@ucl.ac.uk), [John Ashburner](mailto:j.ashburner@ucl.ac.uk)  
MIUA 2019  
https://arxiv.org/abs/1810.03422

- **A Tool for Super-Resolving Multimodal Clinical MRI.**  
[Mikael Brudfors](mailto:brudfors@gmail.com), [Yaël Balbastre](mailto:y.balbastre@ucl.ac.uk), [Parashkev Nachev](mailto:p.nachev@ucl.ac.uk), [John Ashburner](mailto:j.ashburner@ucl.ac.uk)  
Preprint  
https://arxiv.org/abs/1909.01140

The reweighted least squares scheme used in this implementation is 
described in :

- **Joint Total Variation ESTATICS for Robust Multi-Parameter Mapping.**  
[Yaël Balbastre](mailto:y.balbastre@ucl.ac.uk), [Mikael Brudfors](mailto:brudfors@gmail.com), [Michela Azzarito](mailto:michela.azzarito@balgrist.ch ), [Christian Lambert](mailto:christian.lambert@ucl.ac.uk), [Martina F. Callaghan](mailto:m.callaghan@ucl.ac.uk), [John Ashburner](mailto:j.ashburner@ucl.ac.uk)  
Preprint  

## License

This software is released under the 
[GNU General Public License version 3](LICENSE) (GPL v3). As a result, 
you may copy, distribute and modify the software as long as you track 
changes/dates in source files. Any modifications to or software including 
(via compiler) GPL-licensed code must also be made available under the 
GPL along with build & install instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
