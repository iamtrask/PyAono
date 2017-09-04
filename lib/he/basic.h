#include "keys.h"

class ciphertext{
public:
    cipher_text* value;
    int degree;
    public_key* pk;
    
    ciphertext(){};
    
    ~ciphertext(){};
    
    ciphertext(public_key* pk){
        this->degree = 0;
        this->pk = pk;
    }
    
    ciphertext(GEN m, public_key* pk){
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(m);
    }
    
    void initialize(GEN m, public_key* pk){
        this->degree = 2;
        this->pk = pk;
        this->value = pk->encrypt(m);
    }
    
    void initialize(public_key* pk){
        this->degree = 0;
        this->pk = pk;
    }
    
    ciphertext operator+(ciphertext &ct){
        ciphertext result;
        result.value = addition(this->value, ct.value);
        result.pk = this->pk;
        return result;
    }
    
    ciphertext operator*(ciphertext &ct){
        ciphertext result;
        result.value = multiplication(this->value, ct.value);
        result.pk = this->pk;
        return result;
    }
    
    ciphertext operator-(ciphertext &ct){
        ciphertext result;
        result.value = subtraction(this->value, ct.value);
        result.pk = this->pk;
        return result;
    }
    
    GEN decrypt(secret_key sk){
        GEN m = sk.decrypt(this->value);
        return m;
    }
};
