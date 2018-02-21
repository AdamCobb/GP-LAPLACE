# GP-LAPLACE

A Gaussian process based technique for locating attractors from trajectories in time-varying fields.

This repository contains code used in the experiments of our paper: "Identifying Sources and Sinks in the Presence of Multiple Agents with Gaussian Process Vector Calculus" by Adam D. Cobb, Richard Everett, Andrew Markham, and Stephen J. Roberts.

## Abstract 
In systems of multiple agents, identifying the cause of observed agent dynamics is challenging. Often, these agents operate in diverse, non-stationary environments, where models rely on handcrafted environment-specific features to infer influential regions in the system’s surroundings. To overcome the limitations of these inflexible models, we present *GP-LAPLACE*, a technique for locating sources and sinks from trajectories in time-varying fields. Using Gaussian processes, we jointly infer a spatio-temporal vector field, as well as canonical vector calculus operations on that field. Notably, we do this from only agent trajectories without requiring knowledge of the environment, and also obtain a metric for denoting the significance of inferred causal features in the environment by exploiting our probabilistic method. To evaluate our approach, we apply it to both synthetic and real-world GPS data, demonstrating the applicability of our technique in the presence of multiple agents, as well as its superiority over existing methods.

## Example
GP-LAPLACE applied to pelagic seabirds flying over the Mediterranean sea:

![Alt Text](https://github.com/AdamCobb/GP-LAPLACE/blob/master/skl_velocity_Z5.gif)

## Reproducing Results
We have created a number of Jupyter notebooks to reproduce our results:

- `exp_synthetic_stationary.ipynb` - Stationary attractors: Simulating agents in a non-stationary potential field.
- `exp_synthetic_varying-strength.ipynb` - Varying-strength attractors: Simulating agents in a potential field with attractors that vary their strength with time. 
- `exp_synthetic_rotating.ipynb` - Rotating attractors: Simulating agents in a potential field with attractors that change location with time.

Each experiment should take no longer than 20 minutes from start to finish (with load_dict=False). Tested on:
- Ubuntu 16.04, 16GB memory, CPU @ 2.60GHz x 8
- Mac, 8GB memory, CPU @ 3.10GHz x 2

## Getting Started

### Requirements
- [Python==3.5](https://www.python.org/getit/)
- [Tensorflow==1.5](https://www.tensorflow.org/)
- [GPFlow==1.1](https://github.com/GPflow/GPflow)
- [Jupyter](http://jupyter.org)

### Installation
1. Clone GP-LAPLACE and install requirements.
```
cd <installation_path_of_your_choice>
git clone https://github.com/AdamCobb/GP-LAPLACE
cd GP-LAPLACE
pip install -r requirements.txt
```

2. Download and install GPFlow.
```
git clone https://github.com/GPflow/GPflow
pip install GPflow/.
```
3. Run notebooks.
```
cd notebooks
jupyter notebook
```

## Data
> Pollonara, E., Luschi, P., Guilford, T., Wikelski, M., Bonadonna, F. and Gagliardo, A., 2015. Olfaction and topography, but not magnetic cues, control navigation in a pelagic seabird: displacements with shearwaters in the Mediterranean Sea. Scientific reports, 5, p.16486.

## Contact Information
Adam Cobb: acobb@robots.ox.ac.uk

Richard Everett: richard@robots.ox.ac.uk
