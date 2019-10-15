# STB_FVBMM
Video analysis binary features

Instructions
1. Compile the mex files inside LIB/MEX/FV
2. Download the STB features (see DATA folder)
3. Open Main and press F5

This is the code of the paper:

<pre>@ARTICLE{RLeyva19,
author={R. {Leyva} and V. {Sanchez} and C. {Li}},
journal={IEEE Transactions on Image Processing},
title={Compact and Low-Complexity Binary Feature Descriptor and Fisher Vectors for Video Analytics},
year={2019},
volume={28},
number={12},
pages={6169-6184},
abstract={In this paper, we propose a compact and low-complexity binary feature descriptor for video analytics. Our binary descriptor encodes the motion information of a spatio-temporal support region into a low-dimensional binary string. The descriptor is based on a binning strategy and a construction that binarizes separately the horizontal and vertical motion components of the spatio-temporal support region. We pair our descriptor with a novel Fisher Vector (FV) scheme for binary data to project a set of binary features into a fixed length vector in order to evaluate the similarity between feature sets. We test the effectiveness of our binary feature descriptor with FVs for action recognition, which is one of the most challenging tasks in computer vision, as well as gait recognition and animal behavior clustering. Several experiments on the KTH, UCF50, UCF101, CASIA-B, and TIGdog datasets show that the proposed binary feature descriptor outperforms the state-of-the-art feature descriptors in terms of computational time and memory and storage requirements. When paired with FVs, the proposed feature descriptor attains a very competitive performance, outperforming several state-of-the-art feature descriptors and some methods based on convolutional neural networks.},
keywords={computer vision;feature extraction;image motion analysis;image recognition;pattern clustering;vectors;video signal processing;video analytics;binary descriptor;spatio-temporal support region;low-dimensional binary string;vertical motion components;binary data;binary features;binary feature descriptor;Fisher vector scheme;action recognition;computer vision;gait recognition;animal behavior clustering;Feature extraction;Trajectory;Task analysis;Optical imaging;Encoding;Image coding;Computational complexity;Binary features;video analysis;fisher vectors;CNN},
doi={10.1109/TIP.2019.2922826},
ISSN={},
month={Dec},}
</pre>
