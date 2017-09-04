#include "basic.h"

#define n 100
#define lambda 100
#define s 8
#define sigma 4
#define degree_p 7
#define l 64

int test(){
    pari_init(2000000000, 2);
    key_gen keyGen;
    key_pair keys;
    keys = keyGen.generate_key(lambda, l, n, s, sigma, degree_p);
    public_key pk = keys.pk;
    secret_key sk = keys.sk;
    ciphertext ct[4];
    
    GEN m1 = zeromatcopy(1, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=1; j++){
            gel(gel(m1, i), j) = stoi(3);
        }
    }
    GEN m2 = zeromatcopy(1, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=1; j++){
            gel(gel(m2, i), j) = stoi(2);
        }
    }
    GEN m3 = zeromatcopy(1, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=1; j++){
            gel(gel(m3, i), j) = stoi(5);
        }
    }
    
    ct[0].initialize(m1, &pk, pk.params);
    ct[1].initialize(m2, &pk, pk.params);
    ct[2].initialize(m3, &pk, pk.params);
    ct[3] = ((ct[0] + ct[1]) + ct[0]) * ct[0];
    print_GEN(ct[3].decrypt(sk));
    pari_close();
    return 0;
}
