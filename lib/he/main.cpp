#include "basic.h"

#define n 100
#define lambda 100
#define s 8
#define sigma 4
#define degree_p 7
#define l 64

int main(){
    pari_init(2000000000, 2);
    key_gen keyGen;
    key_pair keys;
    keys = keyGen.generate_key(lambda, l, n, s, sigma, degree_p);
    public_key pk = keys.pk;
    secret_key sk = keys.sk;
    ciphertext ct[3];
    ct[0].initialize(zeromatcopy(1, l), &pk);
    ct[1].initialize(zeromatcopy(1, l), &pk);
    ct[2] = ct[0] + ct[1];
    print_GEN(ct[2].decrypt(sk));
    pari_close();
    return 0;
}
