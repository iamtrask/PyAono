#include "keys.h"

class ciphertext{
public:
    cipher_text* value;
    int degree;
    public_key* pk;
    parameters* params;
    
    ciphertext(){};
    
    ~ciphertext(){};
    
    ciphertext(public_key* pk, parameters* params){
        this->degree = 0;
        this->pk = pk;
        this->params = params;
    }
    
    ciphertext(pari_GEN m, public_key* pk, parameters* params){
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(m);
        this->params = params;
    }
    
    ciphertext(int m, public_key* pk, parameters* params){
        pari_GEN inp;
        inp.value = cgetg(2, t_VEC);
        gel(inp.value, 1) = stoi(m);
        pari_GEN messagepacked = create_message_matrix(inp, params->l);
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(messagepacked);
        this->params = params;
    }
    
    ciphertext(cipher_text* ct, public_key* pk, parameters* params){
        this->degree = 2;
        this->pk = pk;
        this->value = ct;
        this->params = params;
    }
    
    void packing_method(pari_GEN m, public_key* pk, parameters* params){
        pari_GEN messagepacked = create_message_matrix(m, params->l);
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(messagepacked);
        this->params = params;
    }
    
    void initialize(pari_GEN m, public_key* pk, parameters* params){
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(m);
        this->params = params;
    }
    
    void initialize(public_key* pk, parameters* params){
        this->degree = 0;
        this->pk = pk;
        this->params = params;
    }
    
    ciphertext operator+(ciphertext &ct){
        ciphertext result;
        result.value = addition(this->value, ct.value, params, pk->pk, pk->g);
        result.params = params;
        result.pk = this->pk;
        return result;
    }
    
    ciphertext operator*(ciphertext &ct){
        ciphertext result;
        result.value = multiplication(this->value, ct.value, params);
        result.params = params;
        result.pk = this->pk;
        return result;
    }
    
    ciphertext operator*(int input){
        ciphertext result;
        pari_GEN inp;
        inp.value = stoi(input);
        result.value = plaintext_multiplication(this->value, inp);
        result.params = params;
        result.pk = this->pk;
        return result;
    }
    
    ciphertext operator-(ciphertext &ct){
        ciphertext result;
        result.value = subtraction(this->value, ct.value, params, pk->pk, pk->g);
        result.params = params;
        result.pk = this->pk;
        return result;
    }
    
    pari_GEN decrypt(secret_key sk){
        pari_GEN m = sk.decrypt(this->value);
        return m;
    }
};

class updation_key{
private:
    pari_GEN X, Y;
    public_key* pk;
    
public:
    parameters* params;
    parameters* params_old;
    globalvars* g;
    
    updation_key(){};
    
    updation_key(pari_GEN X, pari_GEN Y, parameters* params, parameters* params_old, globalvars* g, public_key* pk){
        this->X = X;
        this->Y = Y;
        this->params = params;
        this->params_old = params_old;
        this->g = g;
        this->pk = pk;
    }
    
    void initialize(pari_GEN X, pari_GEN Y, parameters* params, parameters* params_old, globalvars* g, public_key* pk){
        this->X = X;
        this->Y = Y;
        this->params = params;
        this->params_old = params_old;
        this->g = g;
        this->pk = pk;
    }
    
    ciphertext cipher_switch(ciphertext ct){
        GEN A1 = pk->pk->A.value;
        GEN P1 = pk->pk->P.value;
        GEN c1 = ct.value->comp1.value;
        GEN c2 = ct.value->comp2.value;
        int kappa = params->lambda+1;
        GEN f1, f2, f3, E0, F, cdash;
        f1 = zeromatcopy(1, params->n);
        for(int i = 1; i <= params->n; i++){
            for(int j=1; j<=1; j++){
                gel(gel(f1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(params->s)));
            }
        }
        f2 = zeromatcopy(1, params->n);
        for(int i = 1; i <= params->n; i++){
            for(int j=1; j<=1; j++){
                gel(gel(f2, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(params->s)));
            }
        }
        f3 = zeromatcopy(1, params->l);
        for(int i = 1; i <= params->l; i++){
            for(int j=1; j<=1; j++){
                gel(gel(f3, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(params->s)));
            }
        }
        E0 = gadd(gmul(f1, appendmat(A1, P1, params->n, params->l, params->n)), gmul(params->p.value, appendmat(f2, f3, params->n, params->l, 1)));
        GEN bitsc1 = bits(c1, kappa, params_old->n);
        GEN tempcal = gmul(bitsc1, X.value);
        
        F = appendmat(tempcal, gadd(gmul(bits(c1, kappa, params_old->n), Y.value), c2), params->n, params->l, 1);
        
        cdash = gadd(E0, F);
        cipher_text* ctret = new cipher_text;
        pari_GEN cdashret, comp2ret;
        cdashret.value = cdash;
        comp2ret.value = stoi(0);
        ctret->comp1 = cdashret;
        ctret->comp2 = comp2ret;
        ctret->flag = 3;
        ciphertext cp = ciphertext(ctret, pk, params);
        return cp;
    }
    
    void serialize(){
        return;
    }
    
    ~updation_key(){};
};

class updation_key_gen{
public:
    updation_key_gen(){};
    
    updation_key generate_key(key_pair* key1, key_pair* key2){
        
        secret_key sk1 = key1->sk;
        secret_key sk2 = key2->sk;
        public_key pk1 = key1->pk;
        public_key pk2 = key2->pk;
        parameters* params1 = pk1.params;
        parameters* params2 = pk2.params;
        GEN S1, S2, A1, A2, P1, P2;
        int n1, n2, s1, s2;
        S1 = sk1.sk.value;
        S2 = sk2.sk.value;
        A1 = pk1.pk->A.value;
        A2 = pk2.pk->A.value;
        P1 = pk1.pk->P.value;
        P2 = pk2.pk->P.value;
        n1 = params1->n;
        s1 = params1->s;
        n2 = params2->n;
        s2 = params2->s;
        int kappa = params1->lambda+1;
        
        GEN X, E, Y;
        X = zeromatcopy(itos(gmul(stoi(n1),stoi(kappa))), n2);
        long long int nkappa = itos(gmul(stoi(n1),stoi(kappa)));
        for(int i = 1; i <= n2; i++){
            for(int j=1; j<=nkappa; j++){
                gel(gel(X, i), j) = gmodulo(generate_random(params1->lambda).value, params1->q.value);
            }
        }
        E = zeromatcopy(nkappa, params1->l);
        for(int i = 1; i <= params1->l; i++){
            for(int j=1; j<=nkappa; j++){
                gel(gel(E, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params2, pk1.g)), stoi(s2)));
            }
        }
        
        Y = gsub(gadd(gmul(params1->p.value, E), power2(S1, n1, kappa, params1->l, params1->q.value)), gmul(X, S2));
        pari_GEN xret, yret;
        xret.value = X;
        yret.value = Y;
        updation_key updated_key = updation_key(xret, yret, params2, params1, pk1.g, &key2->pk);
        return updated_key;
    }
};

