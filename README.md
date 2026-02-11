# Oddball RNN

Computational modeling and analysis framework for training rate-based recurrent neural networks (RNNs) on oddball / go–nogo style tasks, with downstream population analyses including decoding and lesioning.

Please see https://www.biorxiv.org/content/10.1101/2025.04.09.648012v1 for more details.

## Installation
Clone this repository:
```
git clone https://github.com/NuttidaLab/rnn_oddball.git
cd rnn_oddball
```

Create a conda environment using the tf.yml file included in this repository.
```
conda env create -f tf2.yml
conda activate tf2
```
This environment uses Python 3.6 and TensorFlow 2.6.2. For additional details on training requirements, please refer to the spikeRNN repository: https://github.com/rkim35/spikeRNN

Installation typically takes ~3–4 minutes on a standard workstation (e.g., 4–8 CPU cores, ≥16 GB RAM) but may vary depending on internet speed and package download times.

## Training an RNN on the oddball task
The primary training script is `code/rate/main_sequential.py`

This script implements staged training in which task statistics (e.g., stimulus probability) evolve across training phases, enabling the study of adaptive network dynamics under changing environmental structure.
The sequential training pipeline typically involves:
- Stage 0: High Go probability
- Stage 1: Intermediate probability
- Stage 2: Low Go probability

### Example training script
An example shell launcher is provided:
```
code/rate/go-nogo2.sh
```

This script demonstrates how to configure and train a single model instance, including hyperparameters, task settings, and output directories.
To train a model:
```
cd code/rate
bash go-nogo2.sh
```
NOTE: On a standard workstation setup, training time will vary depending on CPU/GPU availability and network size. On a multi-core AMD CPU with NVIDIA A40 GPU support, training a single RNN with 200 units typically takes ~10–15 minutes to complete across all sequential stages.
