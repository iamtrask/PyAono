# PyAono
A Python implementation of the homomorphic encryption scheme by Yoshinori Aono et al. 

This scheme was first introduced in the paper - "Fast and Secure Linear Regression and Biometric Authentication with Security Update" by Yoshinori Aono et al. (https://pdfs.semanticscholar.org/73f0/aa4e1b47b55f0f3d8464f61750e559067c56.pdf). A really unique feature supported by this proposal is the support for key switching. All the code is written in PARI library (http://pari.math.u-bordeaux.fr/) in C++.

This code also has an independent header containing a PARI implementation of the Knuth-Yao's Algorithm for sampling from a discrete Gaussian distribution. The paper followed for understanding and implementing the algorithm - High Precision Discrete Gaussian Sampling on FPGAs (https://www.esat.kuleuven.be/cosic/publications/article-2372.pdf).

NOTE : Running the homomorphic multiplication gives a ![equation](http://latex.codecogs.com/gif.latex?%24l*l%24) matrix as opposed to the message being ![equation](http://latex.codecogs.com/gif.latex?%241*l%24). This is because after homomorphic multiplication, we don't get message but we get ![equation](http://latex.codecogs.com/gif.latex?%24m%5ET*m%24).

### Importing API
The API can be imported using the command ```import Aono```. It currently supports the following functions and classes:

--------

## Functions Supported:
### 1. ```pari_init(pari_size, max_prime)```
```pari_init()``` is the function that needs to be called before dealing with this API. `pari_size` defines the size of stack we'll be using, and `max_prime` defines the pre computed prime table. 
#### Arguments: ```pari_size (int)```, ```max_prime (int)````
### 2. ```pari_close()```
```pari_close()``` function has to be called at the end of each program to clear the memory used.
### 3. ```create_GEN(x)```
```create_GEN()``` function converts integer `x` to `GEN`
### 4. ```get_element(x, i)```
```get_element()``` function returns the `i^th` element of `GEN` variable `x`
### 5. ```print_GEN(x)```
```print_GEN()``` function prints the `GEN` variable `x`
### 6. ```create_message_matrix(*x / x, l)```
```create_message_matrix(*x / x, l)``` function is used to make a PARI `t_MATRIX` object, which is just a conventional matrix, having a dimension of ```1*l``` (a row matrix) with the elements same as the input vector `x`. If input is just a number, then the first element of the row matrix is set as `x` and rest as `0`.
#### Arguments: ```input (int* or int)```, ```l (int)````
#### 7. ```see_ciphertext(c, index)```
```see_ciphertext(c, index)``` function is used to access the contents of a ciphertext (type `struct cipher_text`). Since the ciphertext has 2 components, `index` is used to decide which component to access, with `0` for component 1 and component 2 otherwise.
#### Arguments: ```c (cipher_text*)```, ```index (int)```

--------
