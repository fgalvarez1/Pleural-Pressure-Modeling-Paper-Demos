# Welcome to the demos for: Álvarez-Barrientos et al. (2026). Physics- and Data-driven reduced-order modeling of spatiotemporal pleural pressure distribution. (Submitted)


<!-- Demos can be browsed statically but also interactively—to start a session and run the code just click on the rocket icon at the top of a tutorial page and then click on "Binder", or directly click on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fgalvarez1/Pleural-Pressure-Modeling-Paper-Demos/master?urlpath=lab/tree/./demos).

## Local run via Docker

Because these demos are computationally intensive, Binder may run out of memory or time out, but you can run them locally using Docker (should work for the foreseeable future):
- Ensure you have [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed and running on your machine.
- Pull the Image:
```bash
docker pull ghcr.io/nchibli/curl-based-reconstruction-paper-demos:latest
```
- Start the Container:
```bash
docker run -p 8888:8888 ghcr.io/nchibli/curl-based-reconstruction-paper-demos:latest
```
- Access the Notebooks: Look at your terminal output for a URL that starts with `http://127.0.0.1:8888/?token=...`, copy and paste that entire link into your web browser to access the demos.

## Local run via system install

You can also run the demos locally by setting up your system (might break at some point):
- Ensure you have [Miniconda](https://docs.anaconda.com/free/miniconda) installed on your system, as well as [git](https://git-scm.com/install).
- Clone the Repository:
```bash
git clone https://github.com/nchibli/Curl-Based-Reconstruction-Paper-Demos.git
cd Curl-Based-Reconstruction-Paper-Demos
```
- Create the conda environment:
```bash
conda env create -f repo2docker/environment.yml
```
This only needs to be done once. Now, every time you want to run the demos, you need to:
- Activate the environment: 
```bash
conda activate notebook
```
- Launch Jupyter:
```bash
jupyter notebook
```
Your default web browser will automatically open.
From there, navigate into the demos/ folder and open the notebooks to get started! -->
