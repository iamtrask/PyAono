# PyAono
A Python implementation of the homomorphic encryption scheme by Yoshinori Aono et al. 

This scheme was first introduced in the paper - <a href="https://pdfs.semanticscholar.org/73f0/aa4e1b47b55f0f3d8464f61750e559067c56.pdf">"Fast and Secure Linear Regression and Biometric Authentication with Security Update"</a> by Yoshinori Aono et al. A really unique feature supported by this proposal is the support for key switching. All the code is written in <a href="http://pari.math.u-bordeaux.fr/">PARI library</a> in C++.

This code also has an independent header containing a PARI implementation of the Knuth-Yao's Algorithm for sampling from a discrete Gaussian distribution. The paper followed for understanding and implementing the algorithm - <a href="https://www.esat.kuleuven.be/cosic/publications/article-2372.pdf">High Precision Discrete Gaussian Sampling on FPGAs</a>.

NOTE : Running the homomorphic multiplication gives a ![equation](http://latex.codecogs.com/gif.latex?%24l*l%24) matrix as opposed to the message being ![equation](http://latex.codecogs.com/gif.latex?%241*l%24). This is because after homomorphic multiplication, we don't get message but we get ![equation](http://latex.codecogs.com/gif.latex?%24m%5ET*m%24).

--------
### What's new in the latest version?
We just rolled out the support for key switching/rotation. Check it out by running `init_key_switching.py`. You can also run `init.py` to check out homomorphic operations working.<br />
This scheme has an `advantage` over other schemes that the key switching can both, increase as well as decrease the security level of the ciphertext! 

### News
New here? Check out our <a href="https://www.youtube.com/watch?v=jKKafEBx5P4">sample tutorial</a> on YouTube recorded by <a href="https://github.com/mayank0403">Mayank</a>. The tutorial covers a basic overview of how to use the API. More in-depth tutorials coming soon. Keep checking out YouTube <a href="https://www.youtube.com/channel/UCzoUqDE_OzYo6lGXtsEbOxQ">channel</a> and our GitHub repositories.

### Contribute
Welcome! If you are intrigued by our work and want to contribute, then please go through the tutorial as mentioned in the `News` section and contact any of the contributors. You should also go through this <a href="https://www.youtube.com/watch?v=sXFmKquiVnk">introductory video</a> by Andrew Trask, of the OpenMined platform.

### Difference between Key Rotation and Security Update
`Key Rotation` means that we change ciphertext associated with the key-pair `(pk1, sk1)`, where `pk1` is used to encrypt the message and `sk1` is used to decrypt to message to another ciphertext associated with keypair `(pk2, sk2)`, without decrypting the first ciphertext.<br />
We specifically call it `Key Rotation` when the parameters have this relation ```n1 = n2```.<br />
We call it `Security Update` if the paramters have this relation ```n1 < n2```, meaning that the security of new ciphertext has been updated and increased.

### Notes regarding the use and order of homomorphic function and key switching [`Important`]
The following cases of operations have been implemented. Support for other cases might require some more literature reading and will follow shortly. Please take a look at he following options which are presently supported: <br />
#### `NOTE`: 
The scheme only supports `ONE` homomorphic multiplication at the moment.<br />
#### `Notations`: 
 * `normal-ciphertext` means a ciphertext which has never been refreshed via either key rotation or security update and has only been operated additively or has never undergone any operation is just fresh.
 * `mult-ciphertext` means a ciphertext which has undergone `ONE` multiplication no matter if it has undergone key rotation or security update before or not. Note that a `mult-ciphertext` `CANNOT` be rotated or refreshed by key rotation or security update because of it's peculiar structure. It has to be decrypted.
 * `rotated-ciphertext` means a ciphertext which has undergone refreshment via either key rotation or security update. <br />
#### `NOTE`: 
If `+` is supported it implies `-` is permitted as well. <br />
#### `Permitted operations and their outputs`:
 * ``` normal-ciphertext + normal-ciphertext = normal-ciphertext```  and ``` normal-ciphertext * normal-ciphertext = mult-ciphertext```
 * ```mult-ciphertext + mult-ciphertext = mult-ciphertext```
 * ```rotated-ciphertext + rotated-ciphertext = rotated-ciphertext``` and ```rotated-ciphertext * rotated-ciphertext = mult-ciphertext```
 * ```rotated-ciphertext + normal-ciphertext = rotated-ciphertext``` and ```rotated-ciphertext * normal-ciphertext = mult-ciphertext```
 * ```normal-ciphertext + rotated-ciphertext = rotated-ciphertext``` and ```normal-ciphertext * rotated-ciphertext = mult-ciphertext```
 * ```mult-ciphertext + normal-ciphertext = mult-ciphertext``` and ```normal-ciphertext + mult-ciphertext = mult-ciphertext```[?]
 * ```mult-ciphertext + rotated-ciphertext = mult-ciphertext``` and ```rotated-ciphertext + mult-ciphertext = mult-ciphertext```[?]

#### `[?] Meaning`:
The above bullets marked with [?] means that these operations don't exist conventionally in the actual paper. I have implemented them as my own decision given they might come in handy during the course of our project. I implemented them in this way (illustrating just one example, rest are similar to this one):<br />
```mult-ciphertext + (normal-ciphertext-of-1s * normal-ciphertext) = mult-ciphertext```<br />
`normal-ciphertext-of-1s` mean that it is a ciphertext which is deliberately made my encrypting a vector of just all 1s, so that if it is multiplied with anything (say `x`), this `x` remains same, but since `x` is of type - `normal-ciphertext` and also this `normal-ciphertext-of-1s` is also of the type - `normal-ciphertext`, this means their homomorphic multiplication `*` will yield a ciphertext encrypting `x` of the type - `mult-ciphertext`. This makes the whole above homomorphic addition as ```mult-ciphertext + mult-ciphertext = mult-ciphertext```.

--------

### Importing API
The API can be imported using the command ```import Aono```. It currently supports the following functions and classes:

--------

## Functions Supported:
### 1. ```pari_init(pari_size, max_prime)```
```pari_init()``` is the function that needs to be called before dealing with this API. `pari_size` defines the size of stack we'll be using, and `max_prime` defines the pre computed prime table. 
#### Arguments: ```pari_size (int)```, ```max_prime (int)```
### 2. ```pari_close()```
```pari_close()``` function has to be called at the end of each program to clear the memory used.

--------

## Classes:
```This class abstracts the GEN variable in C++, making it available through python interface. The class is compatible with +, *, /, -, __getitem__ , %, and print.```<br />

### ```parameters```
  * Class Data:
    1. `n`, `s`, `sigma`, `l`, `lambda` (`ints`)
    2. `q`, `p` (`pari_GEN`)

### ```secret_key```
  * Class Data:
    1. `sk` (`pari_GEN`)
    2. `params` (`parameters*`)

  * ```__init__(self, sk = None, parmas = None)```   
    The constructor initiates class data. 
    Arguments: `sk` (`pari_GEN`), `params` (`parameters*`)

  * ```decrypt(self, ct)```   
    ```decrypt()``` method returns the plaintext (`pari_GEN`) encrypted in ciphertext `ct`.    Arguments: `ct` (`pari_GEN`) 

  * ```serialize(self)```   
    TO BE IMPLEMENTED

### ```public_key```
  * Class Data:
    1. `pk` (`public_key_pack*`)
    2. `params` (`parameters*`)
    3. `g` (`globalvars*`)

  * ```__init__(self, pk = None, params = None, g = None)``` 
    The constructor initiates the class data. 
    Arguments: `pk` (`public_key_pack*`), `params` (`parameters*`), `g` (`globalvars*`)

  * ```encrypt(self, pt)```   
    ```encrypt()``` method returns the ciphertext (`GEN`) which encrypts plaintext `pt`. 
    Arguments: `pt` (`pari_GEN`)

  * ```serialize(self)```   
    TO BE IMPLEMENTED

### ```key_pair```
  * Class Data:
    1. `sk` (`secret_key`)
    2. `pk` (`public_key`)

### ```key_gen```
  * ```generate_key(self, lambda, l, n, s, sigma, degree_p)```
    ```generate_key()``` method returns the keys, which is of type `key_pair`. Here `s` defines the tailprune and `degree_p` defines the `bit_size` of `p` to be generated. 
    Arguments: `lambda` (`int`), `l` (`int`), `n` (`int`), `s` (`int`), `sigma` (`int`), `degree_p` (`int`).

  * ```deserialize(self)```   
    TO BE IMPLEMENTED

### ```ciphertext```   
    The class is compatible with `’+’, '*', and '-' operators`
  * Class Data:
    1. `value` (`pari_GEN`)
    2. `pk` (`public_key*`)
    3. `params` (`parameters*`)

  * ```__init__(self, plaintext = None, pk, params)```  
    The constuctor method takes two arguments: `plaintext` (`pari_GEN` variable), `pk` (`public_key*`), `params` (`parameters*`)

TODO - Remove `params`, it can be taken from `pk`.

  * ```decrypt(self, sk)```  
    ```decrypt()``` method returns the decrypted `ciphertext` which is `pari_GEN` variable. 
Arguments: `sk` (`secret_key*`)

### ```updation_key```
* Class Data:
    1. `pk` (`public_key_pack*`)
    2. `params` (`parameters*`)
    3. `old_params` (`parameters*`)
    4. `g` (`globalvars*`)
    5. `XComponent` (`pari_GEN`)
    6. `YComponent` (`pari_GEN`)

  * ```__init__(self, X, Y, params, params_old, g, pk)``` 
    The constructor initiates the class data. 
    Arguments: `pk` (`public_key_pack*`), `params` (`parameters*`), `g` (`globalvars*`)

  * ```cipher_switch(self, ciphertext)```   
    ```cipher_switch()``` method returns the ciphertext (`ciphertext`) which is obtained by switching the ciphertext passed to function with the updatation key.
    Arguments: `ct` (`ciphertext`)
    
### ```updation_key_gen```

  * ```generate_key(self, key_pair*, key_pair*)```   
    ```generate_key()``` method returns the updation key (`updation_key`) which is the key that can be used to switch any ciphertext which is associated with the first (in order of arguments) `key_pair` to the ciphertext which is associated with the second `key_pair`.
    Arguments: `key1` (`key_pair`), `key2` (`key_pair`)
    
    
--------
   
## TODO:
- [x] Tweaks to homomorphic functions to support rotated or updated ciphertexts. <br/>
- [x] Support for plaintext-ciphertext multiplication. <br/>
- [x] Basic testing suite. <br/>
- [x] Taking care of negative errors by mapping ciphertext to a symmteric group before decryption. <br />
- [x] Add a matrix reading support to testing suite. <br/>
- [x] Add key switching tests to testing suite. <br/>
- [x] Update README with pari_GEN. <br/>

