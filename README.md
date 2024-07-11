# Project-Computational-Neuroscience
Analysis of the Stability of Neuronal Activity in the Frontal Cortex of the Brain 

## 1.Problem statement
 The purpose of this experiment was to determine how short-term memory, which links past emotions to future activities, manifests itself in neuronal dynamics.

## 2.Related works
[1] Inagaki HK, Fontolan L, Romani S, Svoboda K. (2019) Discrete attractor dynamics underlying selective persistent activity in frontal cortex. Nature v566, 212–217. doi: 10.1038/s41586-019-0919-7
[2] Inagaki HK, Inagaki M, Romani S, Svoboda K. (2018) Low-dimensional and monotonic preparatoy activity in mouse anterior lateral motor cortex. Journal of Neuroscience. 3152-17. doi: 10.1523/JNEUROSCI.3152-17.2018

## 3.Proposed Method
In this project, we aim to extract PSTH curves, represent neuron firing in different trials with the same stimulus, and generate tuning curves for neurons in this region. We also seek to compare these curves for different stimuli and functions

## 4.Implementation

### 4.1. Dataset
The data set contains the extracellular and intracellular neurophysiology recordings. These experiments contain data from 31 animals that were trained in a delayed response task. Auditory or whisker tactile cues during the sample epoch instruct the animals to lick one of two possible licking directions. Following a delay epoch, animals report their decisions by directional licking. Duration of the delay epoch was either fixed (fixed delay task: 1.2 s for the intracellular recordings and 2.0 s for the extracellular recordings) or randomized (random delay task) across trials. Recordings were in the anterior motor cortex (ALM) in adult mice. The data represent intracellular and extracellular (64 channels) recordings. For the extracellular recordings, ALM was silenced optogenetically bilaterally during the first half of delay epoch to probe the response of population dynamics.

The data consists of:
Behavior (events related to task and animal responses)
Photostimulation (timing and type)
Intracellular recordings of membrane potential
Processed extracellular recordings (sorted units) 

### 4.2. Model
Gaussian on fixed delay and Fourier on random delay

### 4.3. Evaluate
In the evaluation section, the methods and metrics used to assess the model's performance are detailed
