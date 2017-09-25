#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <pari/pari.h>
#include <time.h>
#define PARI_OLD_NAMES
#include <vector>
#include <string>
#include <sys/time.h>

#define precision 6 // This means 256-bit precision

struct timeval tv;

class pari_GEN{
public:
    GEN value;
    
    pari_GEN(){};
    
    pari_GEN(int x){
        value = stoi(x);
        return;
    }
    
    void initialize(GEN x){
        value = x;
        return;
    }
    
    pari_GEN operator+(const pari_GEN GEN_2){
        pari_GEN result;
        result.value = gadd(this->value, GEN_2.value);
        return result;
    }
    
    pari_GEN operator*(const pari_GEN GEN_2){
        pari_GEN result;
        result.value = gmul(this->value, GEN_2.value);
        return result;
    }
    
    pari_GEN operator/(const pari_GEN GEN_2){
        pari_GEN result;
        result.value = gdiv(this->value, GEN_2.value);
        return result;
    }
    
    pari_GEN operator-(const pari_GEN GEN_2){
        pari_GEN result;
        result.value = gsub(this->value, GEN_2.value);
        return result;
    }
    
    pari_GEN operator%(const pari_GEN GEN_2){
        pari_GEN result;
        result.value = gmodulo(this->value, GEN_2.value);
        return result;
    }
    
    bool operator==(const pari_GEN GEN_2){
        if(!gequal(this->value, GEN_2.value))
            return false;
        else
            return true;
    }
};

struct parameters{
    pari_GEN q, p;
    int lambda;
    int l;
    int s, n; // These are changed during rotation
    int sigma;
};

struct ProbMatrixPack{
    pari_GEN P;
    std::vector<int> startPos;
    bool isInitialized = false;
} *pPackglobal;

struct public_key_pack{
    pari_GEN A;
    pari_GEN P;
    int n, s;
};

struct globalvars{
    ProbMatrixPack* pPack;
    int errorsModulo = 20;
};

struct cipher_text{
    int flag = 1;
    pari_GEN comp1, comp2;
};

struct cipher_text_mult{
    int flag = 2;
    pari_GEN c;
};

GEN get_element(pari_GEN x, int index){
    return gel(x.value, index + 1);
}

GEN get_element(GEN x, int index){
    return gel(x, index + 1);
}

void print_GEN(pari_GEN x){
    printf("%s\n", GENtostr(x.value));
}

void print_GEN(GEN x){
    printf("%s\n", GENtostr(x));
}

GEN create_GEN(int x){
    GEN y = stoi(x);
    return y;
}

GEN create_GEN(char* x){
    GEN y = strtor(x, precision);
    return y;
}

GEN create_GEN(int m, int n){ // Creates a GEN mxn matrix with all 0s
    GEN y = zeromatcopy(m, n);
    return y;
}

double Uniform(void) {
    return ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
}

double Normal(void) {
    return sqrt(-log(Uniform())*2.0) * sin(2.0 * M_PI * Uniform());
}

double Gauss(double mu, double sigma) {
    double z = sqrt(-2.0 * log(Uniform())) * sin(2.0 * M_PI * Uniform());
    return mu + sigma*z;
}

pari_GEN Sample(int n, double sigma) {
    pari_GEN ret;
    ret.value = cgetg(n + 1, t_VEC);
    double z;
    int i;
    
    for (i = 1; i <= n; i++) {
        z = Gauss(0, sigma);
        z = fabs(round(z)); /*absolute value of Gaussian distribution */
        ret.value[i] = (long) stoi((long) z);
    }
    
    return ret;
}

pari_GEN generate_random(int bit_length){
    gettimeofday(&tv, NULL);
    setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
    pari_GEN r;
    r.value = randomi(gshift(gen_1, bit_length));
    return r;
}

GEN getGuassProbability(GEN point, GEN center, parameters* params){
    int sigma = params->sigma;
    GEN twopi = mulir(stoi(2), mppi(precision));
    GEN s = mulir(stoi(sigma), sqrtr(twopi));
    GEN sinv = invr(s);
    if(gcmp(point, strtor("0.00", precision)) == 0)
        return sinv;
    else{
        return gmul(sinv, gexp( gdiv(gneg(gpow(gdiv(gsub(point, center), stoi(sigma)), stoi(2), precision)), strtor("2.00", precision)), precision));
    }
}

ProbMatrixPack* genProbabilityMatrix(parameters* params, char* c){
    int sigma = params->sigma;
    int tailprune = params->s;
    GEN center;
    center = strtor(c, precision);
    
    GEN tempP, beginAddressP, ProbofPoints;
    int bounds = tailprune*sigma;
    ProbofPoints = cgetg(bounds+2, t_REAL);
    int bitprecision = 64*(precision-2);
    tempP = cgetg(bitprecision+1, t_VEC);
    for(int i=1; i<=bitprecision; i++){
        GEN temp = cgetg(bounds+2, t_INT);
        gel(tempP, i) = temp;
    }
    for(int x = bounds; x > 0; x--){
        gel(ProbofPoints, bounds+1-x) = getGuassProbability(gadd(center, stoi(x)), center, params);
    }
    gel(ProbofPoints, bounds+1) = gdiv(getGuassProbability(gadd(center, stoi(0)), center, params), stoi(2));
    
    int i = -1;
    for(int j=0; j<bitprecision; j++){
        GEN temppow = gpow(strtor("2.00", precision), stoi(i), precision);
        i--;
        for(int x = bounds; x >= 0; x--){
            gel(gel(tempP, j+1), bounds+1-x) = stoi(0);
            if(gcmp(gel(ProbofPoints, bounds+1-x), temppow) >= 0){
                gel(gel(tempP, j+1), bounds+1-x) = stoi(1);
                gel(ProbofPoints, bounds+1-x) = gsub(gel(ProbofPoints, bounds+1-x), temppow);
            }
        }
    }
    
    std::vector<int> beginPos;
    
    for(int x = bounds; x >= 0; x--){
        for(int j=0; j<bitprecision; j++){
            if(j == bitprecision-2){
                beginPos.push_back(j);
                break;
            }
            
            if(gcmp(gel(gel(tempP, j+1), bounds+1-x), stoi(1))==0){
                beginPos.push_back(j);
                break;
            }
        }
    }
    ProbMatrixPack* pPack = new ProbMatrixPack;
    pari_GEN P;
    P.value = tempP;
    pPack->P = P;
    pPack->startPos = beginPos;
    pPack->isInitialized = true;
    pPackglobal = pPack;
    return pPack;
    
}

ProbMatrixPack* genProbabilityMatrix(parameters* params, int c){
    int sigma = params->sigma;
    int tailprune = params->s;
    GEN center;
    center = stoi(c);
    
    GEN tempP, beginAddressP, ProbofPoints;
    int bounds = tailprune*sigma;
    ProbofPoints = cgetg(bounds+2, t_REAL);
    int bitprecision = 64*(precision-2);
    tempP = cgetg(bitprecision+1, t_VEC);
    for(int i=1; i<=bitprecision; i++){
        GEN temp = cgetg(bounds+2, t_INT);
        gel(tempP, i) = temp;
    }
    for(int x = bounds; x > 0; x--){
        gel(ProbofPoints, bounds+1-x) = getGuassProbability(gadd(center, stoi(x)), center, params);
    }
    gel(ProbofPoints, bounds+1) = gdiv(getGuassProbability(gadd(center, stoi(0)), center, params), stoi(2));

    
    int i = -1;
    for(int j=0; j<bitprecision; j++){
        GEN temppow = gpow(strtor("2.00", precision), stoi(i), precision);
        i--;
        for(int x = bounds; x >= 0; x--){
            gel(gel(tempP, j+1), bounds+1-x) = stoi(0);
            if(gcmp(gel(ProbofPoints, bounds+1-x), temppow) >= 0){
                gel(gel(tempP, j+1), bounds+1-x) = stoi(1);
                gel(ProbofPoints, bounds+1-x) = gsub(gel(ProbofPoints, bounds+1-x), temppow);
            }
        }
    }
    
    std::vector<int> beginPos;
    
    for(int x = bounds; x >= 0; x--){
        for(int j=0; j<bitprecision; j++){
            if(j == bitprecision-2){
                beginPos.push_back(j);
                break;
            }
            
            if(gcmp(gel(gel(tempP, j+1), bounds+1-x), stoi(1))==0){
                beginPos.push_back(j);
                break;
            }
        }
    }
    ProbMatrixPack* pPack = new ProbMatrixPack;
    pari_GEN P;
    P.value = tempP;
    pPack->P = P;
    pPack->startPos = beginPos;
    pPack->isInitialized = true;
    pPackglobal = pPack;
    return pPack;
    
}

int SampleKnuthYao(int c, parameters* params, globalvars* g){
    int sigma = params->sigma;
    int tailprune = params->s;
    GEN center;
    center = stoi(c);
    
    int bounds, col, d, invsample, pRows, pCols, s, flag, enable, hit;
    unsigned long r;
    bounds = tailprune*sigma;
    d = 0;
    hit = 0;
    invsample = bounds+1;
    
    GEN P = g->pPack->P.value;
    std::vector<int> beginPos = g->pPack->startPos;
    int bitprecision = 64*(precision-2);
    pRows = lg(P)-1;
    pCols = bitprecision;
    
    flag = 1-2*(rand()%2);
    
    int randomBits[pRows];
    int length = sizeof(unsigned long)*8;
    
    
    for(int i=0; i<pRows; i++){
        randomBits[i] = rand()%2;
    }
    
    s = 0;
    enable = 0;
    for(int row = 0; row<pRows; row++){
        if(enable==1)
            break;
        d = 2*d + randomBits[row];
        for(int col = beginPos[row]; col < pCols; col++) {
            d = d - itos(gel(gel(P, col+1), row+1));
            if(d==-1){
                hit = 1;
                s = col;
                enable = 1;
                break;
            }
        }
        
    }
    return (s*flag)+itos(center);
}


globalvars* initialize_sampler(parameters* params){
    globalvars* g = new globalvars;
    g->pPack = genProbabilityMatrix(params, "0.00");
    return g;
}

globalvars* initialize_sampler(parameters* params, int center){
    globalvars* g = new globalvars;
    g->pPack = genProbabilityMatrix(params, center);
    return g;
}

void getGENsize(GEN x){
    printf("%d\n", lg(x)-1);
}

void getGEN_MATRIXsize(GEN x){
    printf("%dx%d\n", lg(gel(x, 1))-1, lg(x)-1);
}

pari_GEN generate_secret_key(parameters* params, globalvars* g){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params);
    }
    GEN q = params->q.value;
    GEN p = params->p.value;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    GEN S;
    S = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(S, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(s)));
        }
    }
    pari_GEN Sret;
    Sret.value = S;
    return Sret;
}

pari_GEN generate_secret_key(parameters* params, int center, globalvars* g){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params, center);
    }
    GEN q = params->q.value;
    GEN p = params->p.value;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    GEN S;
    S = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(S, i), j) = lift(gmodulo(stoi(SampleKnuthYao(center, params, g)), stoi(s)));
        }
    }
    pari_GEN Sret;
    Sret.value = S;
    return Sret;
}

public_key_pack* generate_public_key(pari_GEN sk, parameters* params, globalvars* g){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params);
    }
    GEN q = params->q.value;
    GEN p = params->p.value;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    int modulo = g->errorsModulo;
    
    GEN S = sk.value;
    
    GEN R, A;
    R = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(R, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), stoi(s)));
        }
    }
    A = zeromatcopy(n, n);
    for(int i = 1; i <= n; i++){
        for(int j=1; j<=n; j++){
            gel(gel(A, i), j) = gmodulo(generate_random(params->lambda).value, q);
        }
    }
    
    GEN P, temp;
    temp = RgM_mul(A, S);
    P = gsub(gmul(p, R), temp);

    public_key_pack* pk = new public_key_pack;
    pari_GEN Pret, Aret;
    Pret.value = P;
    Aret.value = A;
    pk->P = Pret;
    pk->A = Aret;
    pk->n = n;
    pk->s = s;
    return pk;
}

public_key_pack* generate_public_key(pari_GEN sk, parameters* params, globalvars* g, int center){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params);
    }
    GEN q = params->q.value;
    GEN p = params->p.value;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    int modulo = g->errorsModulo;
    
    GEN S = sk.value;
    
    GEN R, A;
    R = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(R, i), j) = lift(gmodulo(stoi(SampleKnuthYao(center, params, g)), stoi(s)));
        }
    }
    A = zeromatcopy(n, n);
    for(int i = 1; i <= n; i++){
        for(int j=1; j<=n; j++){
            gel(gel(A, i), j) = gmodulo(generate_random(params->lambda).value, q);
        }
    }
    
    GEN P, temp;
    temp = RgM_mul(A, S);
    P = gsub(gmul(p, R), temp);
    
    public_key_pack* pk = new public_key_pack;
    pari_GEN Pret, Aret;
    Pret.value = P;
    Aret.value = A;
    pk->P = Pret;
    pk->A = Aret;
    pk->n = n;
    pk->s = s;
    return pk;
}


globalvars* set_error_modulo(globalvars* g, int modulo){
    g->errorsModulo = modulo;
    return g;
}

GEN access_value_pk(public_key_pack* pk, char flag){
    if(flag == 'a' || flag == 'A'){
        return pk->A.value;
    }
    else if(flag == 'p' || flag == 'P'){
        return pk->P.value;
    }
    else if(flag == 'n' || flag == 'N'){
        return stoi(pk->n);
    }
    else if(flag == 's' || flag == 'S'){
        return stoi(pk->s);
    }
    else
        return stoi(0);
}

parameters* gen_params(int lambda, int l, int n, int s, int sigma, int degree_p){
    GEN q = nextprime(gpowgs(stoi(2), lambda));
    GEN p = gadd(gpowgs(stoi(2), degree_p), stoi(1));
    parameters* params = new parameters;
    pari_GEN pret, qret;
    pret.value = p;
    qret.value = q;
    params->p = pret;
    params->q = qret;
    params->lambda = lambda;
    params->l = l;
    params->s = s;
    params->n = n; // These are changed during rotation
    params->sigma = sigma;
    return params;
}

pari_GEN create_message_matrix_repeated_input(int , int );
cipher_text* multiplication(cipher_text* , cipher_text* , parameters* );

cipher_text* encrypt_outside_class(pari_GEN m, public_key_pack* pk, parameters* params, globalvars* g){
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

cipher_text* addition(cipher_text* ct_1, cipher_text* ct_2, parameters* params, public_key_pack* pk, globalvars* g){
    if (ct_1->flag == 1 && ct_2->flag == 2){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_1, params);
        GEN ct1, ct2;
        ct1 = gadd(newct->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 2 && ct_2->flag == 1){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_2, params);
        GEN ct1, ct2;
        ct1 = gadd(ct_1->comp1.value, newct->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 2){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_1, params);
        GEN ct1, ct2;
        ct1 = gadd(newct->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 2 && ct_2->flag == 3){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_2, params);
        GEN ct1, ct2;
        ct1 = gadd(ct_1->comp1.value, newct->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 2 && ct_2->flag == 2){
        GEN ct1, ct2;
        ct1 = gadd(ct_1->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 3){
        GEN ct1, ct2;
        ct1 = gadd(ct_1->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 3;
        return ct;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 1){
        GEN q = params->q.value;
        GEN p = params->p.value;
        int lambda = params->lambda;
        int l = params->l;
        int s = params->s;
        int n = params->n;
        int nplusl = n+l;
        GEN c2changed = zeromatcopy(1, nplusl);
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(c2changed, i), j) = gel(gel(ct_2->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(c2changed, i), j) = gel(gel(ct_2->comp2.value, i-n), j);
                }
            }
        }
        
        GEN ct1, ct2;
        ct1 = gadd(ct_1->comp1.value, c2changed);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 3;
        return ct;
    }
    else if (ct_1->flag == 1 && ct_2->flag == 3){
        GEN q = params->q.value;
        GEN p = params->p.value;
        int lambda = params->lambda;
        int l = params->l;
        int s = params->s;
        int n = params->n;
        int nplusl = n+l;
        GEN c1changed = zeromatcopy(1, nplusl);
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(c1changed, i), j) = gel(gel(ct_1->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(c1changed, i), j) = gel(gel(ct_1->comp2.value, i-n), j);
                }
            }
        }
        
        GEN ct1, ct2;
        ct1 = gadd(ct_2->comp1.value, c1changed);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 3;
        return ct;
    }
    else{
        GEN ct1, ct2;
        ct1 = gadd(ct_1->comp1.value, ct_2->comp1.value);
        ct2 = gadd(ct_1->comp2.value, ct_2->comp2.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = ct2;
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        return ct;
    }
}

cipher_text* subtraction(cipher_text* ct_1, cipher_text* ct_2, parameters* params, public_key_pack* pk, globalvars* g){
    if (ct_1->flag == 1 && ct_2->flag == 2){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_1, params);
        GEN ct1, ct2;
        ct1 = gsub(newct->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 2 && ct_2->flag == 1){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_2, params);
        GEN ct1, ct2;
        ct1 = gsub(ct_1->comp1.value, newct->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 2){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_1, params);
        GEN ct1, ct2;
        ct1 = gsub(newct->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 2 && ct_2->flag == 3){
        pari_GEN messagemat = create_message_matrix_repeated_input(1, params->l);
        cipher_text* one_enc = encrypt_outside_class(messagemat, pk, params, g);
        
        cipher_text* newct = multiplication(one_enc, ct_2, params);
        GEN ct1, ct2;
        ct1 = gsub(ct_1->comp1.value, newct->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 2 && ct_2->flag == 2){
        GEN ct1, ct2;
        ct1 = gsub(ct_1->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 2;
        return ct;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 3){
        GEN ct1, ct2;
        ct1 = gsub(ct_1->comp1.value, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 3;
        return ct;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 1){
        GEN q = params->q.value;
        GEN p = params->p.value;
        int lambda = params->lambda;
        int l = params->l;
        int s = params->s;
        int n = params->n;
        int nplusl = n+l;
        GEN c2changed = zeromatcopy(1, nplusl);
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(c2changed, i), j) = gel(gel(ct_2->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(c2changed, i), j) = gel(gel(ct_2->comp2.value, i-n), j);
                }
            }
        }
        
        GEN ct1, ct2;
        ct1 = gsub(ct_1->comp1.value, c2changed);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 3;
        return ct;
    }
    else if (ct_1->flag == 1 && ct_2->flag == 3){
        GEN q = params->q.value;
        GEN p = params->p.value;
        int lambda = params->lambda;
        int l = params->l;
        int s = params->s;
        int n = params->n;
        int nplusl = n+l;
        GEN c1changed = zeromatcopy(1, nplusl);
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(c1changed, i), j) = gel(gel(ct_1->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(c1changed, i), j) = gel(gel(ct_1->comp2.value, i-n), j);
                }
            }
        }
        
        GEN ct1, ct2;
        ct1 = gsub(c1changed, ct_2->comp1.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = stoi(0);
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        ct->flag = 3;
        return ct;
    }
    else{
        GEN ct1, ct2;
        ct1 = gsub(ct_1->comp1.value, ct_2->comp1.value);
        ct2 = gsub(ct_1->comp2.value, ct_2->comp2.value);
        cipher_text* ct = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = ct1;
        comp2ret.value = ct2;
        ct->comp1 = comp1ret;
        ct->comp2 = comp2ret;
        return ct;
    }
}


cipher_text* multiplication(cipher_text* ct_1, cipher_text* ct_2, parameters* params){
    GEN q = params->q.value;
    GEN p = params->p.value;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    int nplusl = n+l;
    
    if (ct_1->flag == 3 && ct_2->flag == 3){
        GEN cbeforemul = ct_1->comp1.value;
        GEN c_1beforemul = ct_2->comp1.value;
        
        GEN cmul = RgM_transmul(cbeforemul, c_1beforemul);
        cipher_text* ret = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = cmul;
        comp2ret.value = stoi(0);
        ret->comp1 = comp1ret;
        ret->comp2 = comp2ret;
        ret->flag = 2;
        return ret;
    }
    else if (ct_1->flag == 3 && ct_2->flag == 1){
        GEN cbeforemul = ct_1->comp1.value;
        GEN c_1beforemul = zeromatcopy(1, nplusl);
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(c_1beforemul, i), j) = gel(gel(ct_2->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(c_1beforemul, i), j) = gel(gel(ct_2->comp2.value, i-n), j);
                }
            }
        }
        
        GEN cmul = RgM_transmul(cbeforemul, c_1beforemul);
        cipher_text* ret = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = cmul;
        comp2ret.value = stoi(0);
        ret->comp1 = comp1ret;
        ret->comp2 = comp2ret;
        ret->flag = 2;
        return ret;
    }
    else if (ct_1->flag == 1 && ct_2->flag == 3){
        GEN cbeforemul = zeromatcopy(1, nplusl);
        GEN c_1beforemul = ct_2->comp1.value;
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(cbeforemul, i), j) = gel(gel(ct_1->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(cbeforemul, i), j) = gel(gel(ct_1->comp2.value, i-n), j);
                    
                }
            }
        }
        
        GEN cmul = RgM_transmul(cbeforemul, c_1beforemul);
        cipher_text* ret = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = cmul;
        comp2ret.value = stoi(0);
        ret->comp1 = comp1ret;
        ret->comp2 = comp2ret;
        ret->flag = 2;
        return ret;
    }
    else{
        GEN cbeforemul = zeromatcopy(1, nplusl);
        GEN c_1beforemul = zeromatcopy(1, nplusl);
        for(int i = 1; i <= nplusl; i++){
            for(int j=1; j<=1; j++){
                if(i<=n){
                    gel(gel(cbeforemul, i), j) = gel(gel(ct_1->comp1.value, i), j);
                    gel(gel(c_1beforemul, i), j) = gel(gel(ct_2->comp1.value, i), j);
                    
                }
                else{
                    gel(gel(cbeforemul, i), j) = gel(gel(ct_1->comp2.value, i-n), j);
                    gel(gel(c_1beforemul, i), j) = gel(gel(ct_2->comp2.value, i-n), j);
                }
            }
        }
        
        GEN cmul = RgM_transmul(cbeforemul, c_1beforemul);
        cipher_text* ret = new cipher_text;
        pari_GEN comp1ret, comp2ret;
        comp1ret.value = cmul;
        comp2ret.value = stoi(0);
        ret->comp1 = comp1ret;
        ret->comp2 = comp2ret;
        ret->flag = 2;
        return ret;
    }
}

pari_GEN create_message_matrix(int message, int l){
    GEN m = zeromatcopy(1, l);
    gel(gel(m, 1), 1) = stoi(message);
    pari_GEN ret;
    ret.value = m;
    return ret;
}

pari_GEN create_message_matrix_repeated_input(int message, int l){
    GEN m = zeromatcopy(1, l);
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=1; j++){
            gel(gel(m, i), j) = stoi(message);
        }
    }
    pari_GEN ret;
    ret.value = m;
    return ret;
}

pari_GEN create_message_matrix(int* message, int l){
    GEN m = zeromatcopy(1, l);
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=1; j++){
            if (i>sizeof(message)/sizeof(*message)){
                gel(gel(m, i), j) = stoi(0);
            }
            else
                gel(gel(m, i), j) = stoi(message[i-1]);
        }
    }
    pari_GEN ret;
    ret.value = m;
    return ret;
}

GEN see_ciphertext(cipher_text* c, int index){
    if(index == 0)
        return c->comp1.value;
    else
        return c->comp2.value;
}

GEN power2(GEN x, int n, int kappa, int l, GEN q){
    GEN power2mat = zeromatcopy(n*kappa, l);
    long long int nkappa = n*kappa;
    GEN pow2 = stoi(1);
    for(int i = 1; i <= l; i++){
        pow2 = stoi(1);
        for(int j=1; j <= kappa; j++){
            for(int k=1; k<=n; k++){
                gel(gel(power2mat, i), (j-1)*n+k) = gmodulo(gmul(gel(gel(x, i), k), pow2), q);
            }
            pow2 = gmul(pow2, stoi(2));
        }
    }
    
    return power2mat;
}

GEN appendmat(GEN m1, GEN m2, int col1, int col2, int row){
    GEN mat = zeromatcopy(row, col1+col2);
    for(int i =1; i<= col1+col2; i++){
        for(int j=1; j<=row; j++){
            if(i<=col1){
                gel(gel(mat, i), j) = gel(gel(m1, i), j);
            }
            else{
                gel(gel(mat, i), j) = gel(gel(m2, i-col1), j);
            }
        }
    }
    
    return mat;
}

GEN bits(GEN m, int kappa, int n){
    GEN mat = zeromatcopy(1, n*kappa);
    long long int nkappa = n*kappa;
    for(int i =1; i<= n; i++){
        GEN bintemp = binary_zv(lift(gel(gel(m, i), 1)));
        bintemp = gtovec(bintemp);
        GEN binx = cgetg(lg(bintemp), t_VEC);
        for(int j=1; j<=lg(bintemp)-1; j++){
            gel(binx, j) = gel(bintemp, lg(bintemp)-j);
        }
        int size = lg(binx)-1;
        for(int j=1; j<=kappa; j++){
            if(j>size){
                gel(gel(mat, i+(j-1)*n), 1) = stoi(0);
            }
            else{
                gel(gel(mat, i+(j-1)*n), 1) = gel(binx, j);
            }
        }
    }
    
    return mat;
}

parameters* get_updation_parameters(parameters* params_old, int n, int s){
    parameters* params = new parameters;
    params->n = n;
    params->s = s;
    params->q = params_old->q;
    params->p = params_old->p;
    params->lambda = params_old->lambda;
    params->l = params_old->l;
    params->sigma = params_old->sigma;
    return params;

}

cipher_text* plaintext_multiplication(cipher_text* ct, pari_GEN input){
    cipher_text* ret = new cipher_text;
    pari_GEN comp1ret, comp2ret;
    comp1ret.value = gmul(ct->comp1.value, input.value);
    comp2ret.value = gmul(ct->comp2.value, input.value);
    ret->comp1 = comp1ret;
    ret->comp2 = comp2ret;
    ret->flag = ct->flag;
    return ret;
}




