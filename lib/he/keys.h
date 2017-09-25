#include "utils.h"

class secret_key{
private:
    pari_GEN sk;
    
public:
    parameters* params;
    
    secret_key(){};
    
    secret_key(pari_GEN sk, parameters* params){
        this->sk = sk;
        this->params = params;
    }
    
    void initialize(pari_GEN sk, parameters* params){
        this->sk = sk;
        this->params = params;
    }
    
    pari_GEN decrypt(cipher_text* ct){
        GEN m;
        if(ct->flag==1)
            m = lift(gmodulo(lift(gadd(gmul(ct->comp1.value, sk.value), ct->comp2.value)), this->params->p.value));
        else if(ct->flag==2){
            GEN SIMatrix = zeromatcopy(params->n+params->l, params->l);
            GEN I = matid(params->l);
            for(int i = 1; i <= params->l; i++){
                for(int j=1; j<=params->l+params->n; j++){
                    if(j<=params->n){
                        gel(gel(SIMatrix, i), j) = gel(gel(sk.value, i), j);
                        
                    }
                    else{
                        gel(gel(SIMatrix, i), j) = gel(gel(I, i), j-params->n);
                    }
                }
            }
            
            m = lift(gmodulo(lift(gmul(RgM_transmul(SIMatrix, ct->comp1.value), SIMatrix)), params->p.value));
            
        }
        else{
            GEN SIMatrix = zeromatcopy(params->n+params->l, params->l);
            GEN I = matid(params->l);
            for(int i = 1; i <= params->l; i++){
                for(int j=1; j<=params->l+params->n; j++){
                    if(j<=params->n){
                        gel(gel(SIMatrix, i), j) = gel(gel(sk.value, i), j);
                        
                    }
                    else{
                        gel(gel(SIMatrix, i), j) = gel(gel(I, i), j-params->n);
                    }
                }
            }
            m = lift(gmodulo(lift(gmul(ct->comp1.value, SIMatrix)), params->p.value));
        }
        pari_GEN mret;
        mret.value = m;
        return mret;
    }
    
    void serialize(){
        return;
    }
    friend class updation_key_gen;
    
    ~secret_key(){};
};

class public_key{
private:
    public_key_pack* pk;
    
public:
    parameters* params;
    globalvars* g;
    
    public_key(){};
    
    public_key(public_key_pack* pk, parameters* params, globalvars* g){
        this->pk = pk;
        this->params = params;
        this->g = g;
    }
    
    void initialize(public_key_pack* pk, parameters* params, globalvars* g){
        this->pk = pk;
        this->params = params;
        this->g = g;
    }
    
    cipher_text* encrypt(pari_GEN m){
        GEN e1 = zeromatcopy(1, params->n);
        GEN e2 = zeromatcopy(1, params->n);
        GEN e3 = zeromatcopy(1, params->l);
        
        for(int i = 1; i <= params->n; i++){
            for(int j=1; j<=1; j++){
                gel(gel(e1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(params->s)));
                gel(gel(e2, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(params->s)));
            }
        }
        for(int i = 1; i <= params->l; i++){
            for(int j=1; j<=1; j++){
                gel(gel(e3, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(params->s)));
            }
        }
        pari_GEN c1, c2;
        c1.value = gadd(RgM_mul(e1, pk->A.value), gmul(params->p.value, e2));
        c2.value = gadd(gadd(RgM_mul(e1, pk->P.value), gmul(params->p.value, e3)), m.value);
        cipher_text* ct = new cipher_text;
        ct->comp1 = c1;
        ct->comp2 = c2;
        return ct;
    }
    
    void serialize(){
        return;
    }

    friend class updation_key_gen;
    friend class updation_key;
    friend class ciphertext;
    
    ~public_key(){};
};

struct key_pair{
    secret_key sk;
    public_key pk;
};

class key_gen{
public:
    key_gen(){};
    
    key_pair generate_key(int lambda, int l, int n, int s, int sigma, int degree_p){
        key_pair keys;
        parameters* params = new parameters;
        pari_GEN p, q;
        q.value = nextprime(gpowgs(stoi(2), lambda));
        p.value = gadd(gpowgs(stoi(2), degree_p), stoi(1));
        params->p = p;
        params->q = q;
        params->lambda = lambda;
        params->l = l;
        params->s = s;
        params->n = n;
        params->sigma = sigma;
        
        globalvars* g = initialize_sampler(params);
        
        pari_GEN temp = generate_secret_key(params, g);
        keys.sk.initialize(temp, params);
        public_key_pack* temp1 = generate_public_key(temp, params, g);
        keys.pk.initialize(temp1, params, g);
        return keys;
    }
    
};

