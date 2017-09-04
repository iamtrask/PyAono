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
    
    ciphertext(GEN m, public_key* pk, parameters* params){
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(m);
        this->params = params;
    }
    
    void initialize(GEN m, public_key* pk, parameters* params){
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
        result.value = addition(this->value, ct.value);
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
    
    ciphertext operator-(ciphertext &ct){
        ciphertext result;
        result.value = subtraction(this->value, ct.value);
        result.params = params;
        result.pk = this->pk;
        return result;
    }
    
    GEN decrypt(secret_key sk){
        GEN m = sk.decrypt(this->value);
        return m;
    }
};
