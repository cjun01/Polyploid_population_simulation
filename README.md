


# Population Simulation Library

A Python library for simulating, visualizing, and analyzing population genomic data.

## install
pip install git+https://github.com/cjun01/population_simulation.git

## Features

- **Hexaploid Population Data Simulation**: Generates artificial data based on given parameters.
- **PCA Visualization**: Visualize the simulated population using PCA.
- **LD Blocks and Intermixing**: Introduces Linkage Disequilibrium blocks and shuffles data between clusters.
- **Linkage Disequilibrium Calculation**: Analyzes the simulated data and plots a heatmap of the Linkage Disequilibrium.

## Installation & Requirements
- **pip install numpy matplotlib scikit-learn pandas seaborn tqdm numba

## How to use
- **Simulating and Visualizing the Population**

from your_module_name import simulate_and_visualize_population

df, data = simulate_and_visualize_population(2000, 15000, 6, 5)

<sub>note:population size = 2000<sub>

<sub>number of marker = 15000<sub>

<sub>ploidy level =6<sub>

<sub>number of population structure =5<sub>

- **Calculating and Visualizing Linkage Disequilibrium**

from your_module_name import calculate_r2

ld_r2 = calculate_r2(data)

note: 

<sub>data is generated from simulate_and_visualize_population step<sub>


## License

This project is open-source and available under the MIT License. See the LICENSE file for more information.
