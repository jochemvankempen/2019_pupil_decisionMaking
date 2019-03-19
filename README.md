These scripts reproduce the analysis in the paper: van Kempen et al.,
(2019) 'Behavioural and neural signatures of perceptual decision-making are modulated by pupil-linked arousal'
DOI: 10.7554/eLife.42541
https://elifesciences.org/articles/42541

These scripts are based on the original scripts for the paper Newman et al. (2017), 
Journal of Neuroscience (http://www.jneurosci.org/content/37/12/3378). https://github.com/gerontium/big_dots

All data can be found at https://figshare.com/s/8d6f461834c47180a444

-------------------------------------------------------------------------

Initial preprocessing analysis steps (e.g. that determine which channels are noisy and need to be interpolated), can be found at https://github.com/gerontium/big_dots. 

The place to start is batch_bigDots.m, from there all other functions are called.
The scripts in this repository reproduce the preprocessing steps such as interpolating, filtering and epoching, and further more detailed analyses described in the paper.

