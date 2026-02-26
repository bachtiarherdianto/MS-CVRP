# Guided Multiple Search algorithm for solving the CVRP

## Overview

This repository contains the implementation of the **Guided Multiple Search (Guided-MS)** metaheuristic algorithm presented in our paper:

[***Feature-Guided Metaheuristic with Diversity Management for Solving the Capacitated Vehicle Routing Problem***](https://doi.org/10.1016/j.ejor.2025.12.029)

For details regarding the data analysis and the "learning from past solutions" framework described in the paper, please refer to the [MS-Feature](https://github.com/bachtiarherdianto/MS-Feature) repository.

## Installation

### Clone the Repository

Start by cloning the repository to your local machine:

```bash
git clone https://github.com/bachtiarherdianto/MS-CVRP.git
```

### Build the Algorithm

Navigate to the project directory and create a build folder:

```bash
cd MS-CVRP
mkdir build && cd build
```

#### Build Options

The build system uses `cmake` and supports several compile-time options:

| Option              | Description                                                            |
| - | - |
| `CMAKE_BUILD_TYPE`  | Set to `Release` for optimized builds.                                 |
| `ENABLE_VERBOSE`    | Prints informational messages during processing (`0` = off, `1` = on). |
| `ENABLE_LOG_OUTPUT` | Logs every new solution found (can slow down execution).               |
| `ENABLE_GUIDANCE`   | Activates the proposed feature-based guidance (`0` = off, `1` = on).   |

#### Build Commands

**1. Baseline MS Algorithm (without guidance):**

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_VERBOSE=1 -DENABLE_LOG_OUTPUT=0 -DENABLE_GUIDANCE=0
make -j
```

**2. With Proposed Feature-Based Guidance:**

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_VERBOSE=1 -DENABLE_LOG_OUTPUT=0 -DENABLE_GUIDANCE=1
make -j
```

> You can adjust these options depending on whether you want verbose output, logs, or the guidance feature enabled.

### Running the Code

After building, you can run the executable by providing a CVRP instance file as the first argument.

#### Example

Assuming the executable is in the `build` folder:

```bash
cd /home/user/git/MS-CVRP/build
./MS-CVRP <path-to-instance> [optional-args]
```

To see all available optional command-line arguments:

```bash
./MS-CVRP --help
```

## Citation

If you find this repository helpful for your research, please cite:

```bibtex
@article{HERDIANTO2025,
    title = {Feature-guided metaheuristic with diversity management for solving the capacitated vehicle routing problem},
    journal = {European Journal of Operational Research},
    year = {2025},
    issn = {0377-2217},
    doi = {https://doi.org/10.1016/j.ejor.2025.12.029},
    url = {https://www.sciencedirect.com/science/article/pii/S0377221725009968},
    author = {Bachtiar Herdianto and Romain Billot and Flavien Lucas and Marc Sevaux},
}
```

## License

This code is provided **for non-commercial research purposes only**. For commercial or business usage, please contact the authors for permission.

## Contact

For questions, collaborations, or updates regarding this repository: [bachtiarherdianto@gmail.com](mailto:bachtiarherdianto@gmail.com)

## Acknowledgment

This research was supported by the French Agence Nationale de la Recherche (ANR), under grant [ANR-20-THIA-0019 (AI@IMT)](https://anr.fr/Projet-ANR-20-THIA-0019) and [ANR-22-CE22-0016-01 (MAMUT)](https://anr.fr/Projet-ANR-22-CE22-0016). We would also like to acknowledge the open-source projects and communities whose tools and contributions made this research possible.

