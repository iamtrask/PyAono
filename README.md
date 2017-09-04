# PyAono
A Python implementation of the homomorphic encryption scheme by Yoshinori Aono et al. 

This scheme was first introduced in the paper - "Fast and Secure Linear Regression and Biometric Authentication with Security Update" by Yoshinori Aono et al. (https://pdfs.semanticscholar.org/73f0/aa4e1b47b55f0f3d8464f61750e559067c56.pdf). A really unique feature supported by this proposal is the support for key switching. All the code is written in PARI library (http://pari.math.u-bordeaux.fr/) in C++.

This code also has an independent header containing a PARI implementation of the Knuth-Yao's Algorithm for sampling from a discrete Gaussian distribution. The paper followed for understanding and implementing the algorithm - High Precision Discrete Gaussian Sampling on FPGAs (https://www.esat.kuleuven.be/cosic/publications/article-2372.pdf).

NOTE : Running the homomorphic multiplication gives a ![equation](http://latex.codecogs.com/gif.latex?%24l*l%24) matrix as opposed to the message being ![equation](http://latex.codecogs.com/gif.latex?%241*l%24). This is because after homomorphic multiplication, we don't get message but we get ![equation](http://latex.codecogs.com/gif.latex?%24m%5ET*m%24).
