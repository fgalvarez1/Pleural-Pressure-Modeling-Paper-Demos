# Demos for: Physics- and Data-driven reduced-order modeling of spatiotemporal pleural pressure distribution

**Reference:** Álvarez-Barrientos et al. (2026). *Physics- and Data-driven reduced-order modeling of spatiotemporal pleural pressure distribution.* (Submitted).

Demos can be browsed statically but also interactively—to start a session and run the code just click on the rocket icon at the top of a tutorial page and then click on "Binder", or directly click on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fgalvarez1/Pleural-Pressure-Modeling-Paper-Demos/HEAD?urlpath=lab/tree/demos).


## Overview of Demos

Depending on the computational resources required, we recommend different methods for running the notebooks:

**1. Lightweight Demos (Run in Binder)**  
These notebooks can be run interactively right in your browser without any local setup.
* [Construction of the pleural pressure reduced-order model](https://fgalvarez1.github.io/Pleural-Pressure-Modeling-Paper-Demos/demos/Fig3-9_pleural_pressure_reduced_model.html)
* [Synthetic examples](https://fgalvarez1.github.io/Pleural-Pressure-Modeling-Paper-Demos/demos/FigAppB_synthetic_example.html)

**2. Computationally Intensive Demos (Run Locally)**  
The following notebooks require significant memory and processing power. Because Binder may run out of memory or time out, we highly recommend running them locally via **Docker** or a **System Install**:
* [Estimating the pleural pressure](https://fgalvarez1.github.io/Pleural-Pressure-Modeling-Paper-Demos/demos/Fig2a_estimate_pleural_pressure.html)
* [Computing transport plans](https://fgalvarez1.github.io/Pleural-Pressure-Modeling-Paper-Demos/demos/Fig2b_compute_transport_plans.html)
* [Saving the results in the target mesh](https://fgalvarez1.github.io/Pleural-Pressure-Modeling-Paper-Demos/demos/Fig2c_save_transported_fields.html)


## Running Locally


### Option A: Docker (Recommended)
The Docker image should work for the foreseeable future.

1. Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is installed and running on your machine.
2. Pull the latest image:
```bash
docker pull ghcr.io/fgalvarez1/pleural-pressure-modeling-paper-demos:latest
```
3. Start the container:
```bash
docker run -p 8888:8888 ghcr.io/fgalvarez1/pleural-pressure-modeling-paper-demos:latest
```
4. **Access the Notebooks:** Look at your terminal output for a URL that starts with `http://127.0.0.1:8888/?token=...`. Copy and paste that entire link into your web browser to access the demos.


### Option B: System Install (Conda)

You can also run the demos locally by setting up the Python environments directly on your system.

1. Ensure you have [Miniconda](https://docs.anaconda.com/free/miniconda) and [Git](https://git-scm.com/install) installed.
2. Clone the repository:
```bash
git clone https://github.com/fgalvarez1/Pleural-Pressure-Modeling-Paper-Demos.git
cd Pleural-Pressure-Modeling-Paper-Demos
```
3. Create the necessary conda environments (this only needs to be done once):
```bash
conda env create -f .repo2docker/environment.yml
conda env create -f .repo2docker/environment_pot.yml
```
4. Every time you want to run the demos, you need to activate the environment. Note that each demo indicates the environment needed (`notebook` or `pot`):
```bash
conda activate notebook
```
or
```bash
conda activate pot
```
Then, launch Jupyter:
```bash
jupyter notebook
```
Your default web browser will automatically open. From there, navigate into the `demos/` folder and open the notebooks to get started!


## Acknowledgements

The Jupyter books are created based on the [Jupyter Book Project](https://github.com/jupyter-book/jupyter-book).