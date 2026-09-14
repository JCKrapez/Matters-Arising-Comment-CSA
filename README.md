# Matters Arising Comment CSA

## Overview

This repository provides MATLAB scripts developed or modified as part of the evaluation work for the paper  
He, Y., Chen, M. K., Huang, M. et al. Dispersive Meta-lens Thermometry for High-temperature Measurements. Nature Communications 16, 10090 (2025). https://doi.org/10.1038/s41467-025-65171-7
and the related MATLAB scripts available in the repository https://github.com/KirinShi/Temperature-reconstruction-using-CSA
The original repository provides MATLAB scripts and example datasets for temperature reconstruction from spectral measurement data according to the CSA method (Chameleon Swarm Algorithm).
A document 'Matters Arising' was submitted to Nature Communications which provides results obtained with the codes in the present repository. All scripts except one come from the original repository https://github.com/KirinShi/Temperature-reconstruction-using-CSA. Some of them were modified in order to :
-	fix the bug detected in CSA ; CSA can also run as originally, with the bug present;
-	display intermediate results (raw radiance signals, their ratio, the inferred emissivity, the emissivity and radiance temperatures obtained after the last iteration, the final mean temperature and the associated emissivity spectrum);
-	display additional information such as the diagram of permitted solutions, bias and RMS errors ;
-	perform data processing as if the material were grey (synthetic radiance signal).

When a script was modified, the suffix '_VMA' (Version Matters Arising) was added to its name.  
The experimental data are those available in the original repository (tests performed on a blackbody and an alumina plate).

In addition, a script was specifically developed: MWT_permitted_solutions_VMA.m. It was used to obtain the results presented in the main paper 'Matters Arising: Comment on Dispersive Meta-lens Thermometry for High-temperature Measurements'. The results from all other scripts are in the Supplementary material. ---

The experimental data obtained by He et al. have to be downloaded from the repository https://github.com/KirinShi/Temperature-reconstruction-using-CSA.

## Repository Structure

```text
Project/
├── Permitted_solutions/
│    └── MWT_permitted_solutions_VMA.m
├── Temperature-reconstruction-using-CSA-VMA/
    ├── src/ 
        ├── process/
        │    ├── temperature_reconstruction_blackbody_cali_VMA.m
        │    └── temperature_reconstruction_materials_VMA.m
        └── temperature_retrieval/
            ├── spectrumTemperature5_5_VMA.m
            ├── temrecon_singlepoint_VMA.m
            └── CSA/
                ├── Chameleon_VMA.m
                ├── get_orthonormal.m
                ├── initialization.m
                ├── rotation.m
                ├── RotMatrix.m

```

### Source Code

> `Temperature-reconstruction-using-CSA-VMA/src/temperature_retrieval/`
> This folder contains the core temperature reconstruction algorithms.
* `CSA/`: Implementation of the CSA optimization solver used in the temperature reconstruction process.
* `Chameleon_VMA.m`: the original source code Chameleon.m contains a bug. The present version allows performing the calculations with the bug fixed (or the bug present, to allow compararison).
* `spectrumTemperature5_5_VMA.m`: Main entry point for temperature reconstruction. It defines the mathematical formulation of the temperature reconstruction problem and performs the reconstruction procedure.
* `temrecon_singlepoint_VMA.m`: MATLAB script for single-point temperature reconstruction.

> `Temperature-reconstruction-using-CSA-VMA/src/process/`
> This folder contains modified processing scripts for the provided example datasets. 
* `temperature_reconstruction_blackbody_cali_VMA.m`: Modified processing script for the blackbody temperature reconstruction datasets. In addition to the original version, many options are added to display intermediate and additional results, with the bug present or the bug fixed.
* `temperature_reconstruction_materials_VMA.m`: Modified processing script for the material temperature reconstruction datasets. In addition to the original version, many options are added to display intermediate and additional results, with the bug present or the bug fixed.

> `Permitted_solutions/`
> This folder supplement the original repository with the following script:
* `MWT_permitted_solutions_VMA.m`: Script to display the permitted temperature and emissivity solutions for two generic cases: blackbody and alumina with emissivities assumed linear or taken from the literature ('Matters Arising' core paper).

### Data (to be downloaded from the repository https://github.com/KirinShi/Temperature-reconstruction-using-CSA)

> `data/` 
* `blackbody/`: Dataset corresponding to the blackbody reconstruction case. The folders 1 to 21 correspond to different temperature conditions from 1873 K to 1673 K with an interval of 10 K.
* `materials/`: Dataset corresponding to the material reconstruction case.

---

## Requirements

* MATLAB R2021a
The original repository requires Matlab R2024a. A few changes were made to require only R2021a for the present repository
No additional MATLAB toolboxes are required unless otherwise specified.

---

## Usage

1. Open MATLAB R2021a.
2. Set the repository root directory as the MATLAB current working directory.
3. Download the folder `data/` from the original repository https://github.com/KirinShi/Temperature-reconstruction-using-CSA and store it in the folder Temperature-reconstruction-using-CSA-VMA/
3. Run the required script:

For displaying the permitted temperature and emissivity solutions for the two generic cases (blackbody and alumina) presented in the core paper "Matters Arising" :

`MWT_permitted_solutions_VMA.m` in `Permitted_solutions/`
To run on your data, replace the values of the emissivity set v_emiss (the temperature value is asked anyway to the user and stored in the variable Temper).

For blackbody temperature reconstruction (as presented in the Supplementary material of the "Matters Arising" document):

`temperature_reconstruction_blackbody_cali_VMA.m` in `Temperature-reconstruction-using-CSA-VMA/src/process/`

For material temperature reconstruction (as presented in the Supplementary material of the "Matters Arising" document):

`temperature_reconstruction_materials_VMA.m` in `Temperature-reconstruction-using-CSA-VMA/src/process/`

The last two scripts will automatically load the corresponding datasets from the `Temperature-reconstruction-using-CSA-VMA/data/` directory and perform the temperature reconstruction.
Executing these scripts is an interactive process, various options are offered for analyzing and displaying intermediate and final results (e.g., with or without the bug correction, for one or more pixels or the entire image, etc.). When the user is prompted to enter a parameter value, information is provided regarding the value used to obtain the results presented in the "Matters Arising" document.
Typical execution time can be up to 7 seconds per pixel.

---

## Third-Party Code

This repository contains third-party codes (CSA and the original scripts and data from https://github.com/KirinShi/Temperature-reconstruction-using-CSA). The original copyright and license information are preserved in the corresponding third-party directories.

---

## License

The original code and datasets provided in this repository are released under the MIT License unless otherwise stated.

Third-party components are distributed under their own licenses.

---

## Citation

If you use this repository in your research, please cite the associated publications:

```text
Krapez, J.-C. Matters Arising: Comment on Dispersive Meta-lens Thermometry for High-temperature Measurements. submitted to Nature Communications (2026).

He, Y., Chen, M. K., Huang, M. et al. Dispersive Meta-lens Thermometry for High-temperature Measurements. Nature Communications 16, 10090 (2025). https://doi.org/10.1038/s41467-025-65171-7

Braik, M. S. Chameleon Swarm Algorithm: A bio-inspired optimizer for solving engineering design problems. Expert Systems with Applications, 174, 114685 (2021). https://doi.org/10.1016/j.eswa.2021.114685 
```

---

## Contact

For questions regarding the modified code, please contact J.-C. Krapez, jean-claude.krapez@onera.fr

