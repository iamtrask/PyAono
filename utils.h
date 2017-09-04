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

#define precision 6 // This means 256-bit precision

struct parameters{
    GEN q, p;
    int lambda;
    int l;
    int s, n; // These are changed during rotation
    int sigma;
};

struct ProbMatrixPack{
    GEN P;
    vector<int> startPos;
    bool isInitialized = false;
} *pPackglobal;

struct public_key{
    GEN A;
    GEN P;
    int n, s;
};

struct globalvars{
    ProbMatrixPack* pPack;
    int errorsModulo = 20;
};

GEN get_element(GEN x, int index){
    return gel(x, index + 1);
}

void print_GEN(GEN x){
    printf("%s\n", GENtostr(x));
}

GEN create_GEN(int x){
    GEN y = stoi(x);
    return y;
}

GEN create_GEN(float x){
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

GEN Sample(int n, double sigma) {
    GEN ret = cgetg(n + 1, t_VEC);
    double z;
    int i;
    
    for (i = 1; i <= n; i++) {
        z = Gauss(0, sigma);
        z = fabs(round(z)); /*absolute value of Gaussian distribution */
        ret[i] = (long) stoi((long) z);
    }
    
    return ret;
}

GEN getGuassProbability(GEN point, GEN center, parameters* params){
    int sigma = params->sigma;
    GEN twopi = mulir(stoi(2), mppi(precision));
    GEN s = mulir(sigma, sqrtr(twopi));
    GEN sinv = invr(s);
    if(gcmp(point, strtor("0.00", precision)) == 0)
        return sinv;
    else{
        return gmul(sinv, gexp( gdiv(gneg(gpow(gdiv(gsub(point, center), sigma), stoi(2), precision)), strtor("2.00", precision)), precision));
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
    
    vector<int> beginPos;
    
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
    pPack->P = tempP;
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
    
    vector<int> beginPos;
    
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
    pPack->P = tempP;
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
    
    GEN P = g->pPack->P;
    vector<int> beginPos = g->pPack->startPos;
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
    g->pPack = genProbabilityMatrix(parameters* params, "0.00");
}

globalvars* initialize_sampler(parameters* params, int center){
    globalvars* g = new globalvars;
    g->pPack = genProbabilityMatrix(parameters* params, center);
}

void getGENsize(GEN x){
    printf("%d\n", lg(x)-1);
}

void getGEN_MATRIXsize(GEN x){
    printf("%dx%d\n", lg(gel(x, 1))-1, lg(x)-1);
}

GEN generate_secret_key(parameters* params, globalvars* g){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params);
    }
    GEN q = params->q;
    GEN p = params->p;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    GEN S;
    S = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(S, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), s));
        }
    }
    return S;
}

GEN generate_secret_key(parameters* params, int center, globalvars* g){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params, center);
    }
    GEN q = params->q;
    GEN p = params->p;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    GEN S;
    S = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(S, i), j) = lift(gmodulo(stoi(SampleKnuthYao(center, params, g)), s));
        }
    }
    return S;
}

public_key* generate_public_key(GEN sk, parameters* params, globalvars* g){
    if(g->pPack->isInitialized == false){
        g = initialize_sampler(params);
    }
    GEN q = params->q;
    GEN p = params->p;
    int lambda = params->lambda;
    int l = params->l;
    int s = params->s;
    int n = params->n;
    int modulo = g->errorsModulo;
    
    GEN S = sk;
    
    GEN R, A;
    R = zeromatcopy(n, l);
    
    for(int i = 1; i <= l; i++){
        for(int j=1; j<=n; j++){
            gel(gel(R, i), j) = lift(gmodulo(stoi(SampleKnuthYao(0, params, g)), s));
        }
    }
    A = zeromatcopy(n, n);
    for(int i = 1; i <= n; i++){
        for(int j=1; j<=n; j++){
            gel(gel(A, i), j) = gmodulo(stoi(rand()%modulo), q);
        }
    }
    
    GEN P, temp;
    temp = RgM_mul(A, S);
    P = gsub(gmul(p, R), temp);

    public_key* pk = new public_key;
    pk->P = P;
    pk->A = A;
    pk->n = n;
    pk->s = s;
    return pk;
}


globalvars* set_error_modulo(globalvars* g, int modulo){
    g->errorsModulo = modulo;
    return g;
}




