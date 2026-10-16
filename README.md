# Code and Data for Formalizing Poisson-Boltzmann Theory for Field-Tunable Nanofluidic Devices

This repository contains the code and data supporting the manuscript:

**Formalizing Poisson-Boltzmann Theory for Field-Tunable Nanofluidic Devices**

Zhongyuan Zhao, Chudi Qi, Yuheng Li, Shoushan Fan, Qunqing Li, and Yang Wei

Department of Physics and Tsinghua-Foxconn Nanotechnology Research Center,  
State Key Laboratory of Low-Dimensional Quantum Physics,  
Tsinghua University, Beijing 100084, China.

preprint at https://doi.org/10.48550/ARXIV.2604.14777 (2026). 

## Overview

This work reformulates the Poisson-Boltzmann theory for ion transport in nanoconfined systems. By introducing two dimensionless parameters describing confinement and surface electric-field strength, the framework identifies distinct electric-double-layer regimes and connects them to experimentally relevant nanofluidic transport behavior.

The framework is applied to field-tunable nanofluidic devices, particularly ionic transistors. It reproduces characteristic conductance-concentration relations, explains gate-controlled transport with different polarities, and identifies a thermodynamic subthreshold-swing limit of 60 mV/dec at room temperature.

## Repository Scope

The materials in this repository support the numerical calculations, analyses, and figures reported in the manuscript. 

## Scientific Framework

The model considers a nanochannel in equilibrium with an electrolyte reservoir. The ionic distributions are described by the Poisson-Boltzmann equation. After nondimensionalization, the system is characterized primarily by:

- **Confinement parameter** `γ`, which describes the channel size relative to the Debye screening length;
- **Relative field-strength parameter** `χ`, which describes the surface electric field relative to ionic screening.

The resulting parameter space separates the confined electric double layer into linear-response, electric-double-layer-overlap, and surface-accumulation regimes. These regimes provide the basis for interpreting and optimizing field-tunable ion transport.

## About the code:

The code.py file contains the code solving the Poisson-Boltzmann equation in a nano-slit configuration. The input variables are the channel height, the Debye length and the Gouy-Chapman length. The solver utilizes an self-adaptive algorithm to solve the PB boundary problem. The solver yields the data of z-coordinates, potential profile and electric field profile, which can be further processed to derive Phi_s, Phi_D, Gamma and Lambda to analyze the field-tunable behavior of a nanofluidic device.

## About the data:

The data supporting the findings in the main text of the manuscript are provided in the data.csv file, which are the mapping data of Phi_s, Phi_D, Gamma and Lambda over the gamma-chi parameter space. They constitute the least necessary set of data to fully reproduce the main results shown in the manuscript. Further information can be provided upon reasonable request.
