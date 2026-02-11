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
