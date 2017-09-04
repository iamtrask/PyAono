//
//  knuthYaoSampler.h
//  AonoEnc
//
//  Created by Mayank Rathee on 25/08/17.
//  Copyright © 2017 Mayank Rathee. All rights reserved.
//

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <pari/pari.h>
#include <time.h>
#include <string.h>
#define PARI_OLD_NAMES
#include <iostream>
#include <vector>

using namespace std;

struct ProbMatrixPack{
    GEN P;
    vector<int> startPos;
} *pPacklocal;

ProbMatrixPack* pPack = new ProbMatrixPack;


//mulir(sinv, gexp( gdiv(gneg(gpow(gdiv(gsub(point, center), sigma), stoi(2), precision)), strtor("2.00", precision)), precision));

GEN getGuassProbability(GEN point, GEN sigma, GEN center, int precision){
    GEN twopi = mulir(stoi(2), mppi(precision));
    GEN s = mulir(sigma, sqrtr(twopi));
    GEN sinv = invr(s);
    
    // if point is 0
    if(gcmp(point, strtor("0.00", precision)) == 0)
        return sinv;
    else{
        return gmul(sinv, gexp( gdiv(gneg(gpow(gdiv(gsub(point, center), sigma), stoi(2), precision)), strtor("2.00", precision)), precision));
    }
}

// Function to generate the probability matrix which is a substitute of the DDT Tree. Some really fine optimizations are not considered here which help prune this matrix

ProbMatrixPack* getProbabilityMatrix(int sigma, char* c, int precision, int tailprune){
    GEN center;
    center = strtor(c, precision);
    
    GEN tempP, beginAddressP, ProbofPoints;
    int bounds = tailprune*sigma;
    // Now Points will run in [-bound, bound]. Since, all the points here except 0 will have 2 occurences, one of positive side and one on negative side (like x and -x). So, the probablility of 0 will be reduced to half later in this function
    // We consider only the positive points along with 0 (with half probability)
    ProbofPoints = cgetg(bounds+2, t_REAL); // extra 1 for 0
    // For a 64-bit machine the number of bits in long is 64
    int bitprecision = 64*(precision-2);
    cout<<"Knuth Yao Probability matrix is moving forward with a bit precision of "<<bitprecision<<endl;
    
    tempP = cgetg(bitprecision+1, t_VEC);
    for(int i=1; i<=bitprecision; i++){
        GEN temp = cgetg(bounds+2, t_INT);
        gel(tempP, i) = temp;
    }
    // Some useless statements
    //cout<<GENtostr(gexp(stoi(-1), precision));
    //cout<<GENtostr(gdiv(center, stoi(2)));
    //cout<<GENtostr(gpow(strtor("2.00", precision), stoi(2), precision));
    //cout<<GENtostr(getGuassProbability(strtor("6.00", precision), stoi(sigma), center, precision));
    
    for(int x = bounds; x > 0; x--){
        gel(ProbofPoints, bounds+1-x) = getGuassProbability(gadd(center, stoi(x)), stoi(sigma), center, precision);
    }
    // We have found the probability for the rest of points except the point x = 0, which will have 1/2 probability because it does not have a negative counterpart.
    gel(ProbofPoints, bounds+1) = gdiv(getGuassProbability(gadd(center, stoi(0)), stoi(sigma), center, precision), stoi(2));
    
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
    
    //cout<<GENtostr(gel(gel(tempP, 1), 9))<<endl;
    
    for(int x = bounds; x >= 0; x--){
        for(int j=0; j<bitprecision; j++){
            //cout<<GENtostr(gel(gel(tempP, j+1), bounds+1-x))<<"\n";
        }
        //cout<<endl<<endl;
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
    
    pPack->P = tempP;
    pPack->startPos = beginPos;
    
    
    
    //pPack = pPacklocal;
    return pPack;
    // Some part remaining
    
}

int SampleKnuthYao(int tailprune, int sigma, int c, int precision){
    GEN center;
    center = stoi(c);
    
    int bounds, col, d, invsample, pRows, pCols, s, flag, enable, hit;
    unsigned long r;
    bounds = tailprune*sigma;
    d = 0;
    hit = 0;
    invsample = bounds+1;
    
    //ProbMatrixPack* pPack = getProbabilityMatrix(4, "3.455555554534535353253425234543534535345245235312345678901234567890", 6, 4);
    
    GEN P = pPack->P;
    vector<int> beginPos = pPack->startPos;
    int bitprecision = 64*(precision-2);
    pRows = lg(P)-1;
    pCols = bitprecision;
    
    flag = 1-2*(rand()%2); // Requires change in the PRNG. A stronger one is required.
    
    int randomBits[pRows];
    int length = sizeof(unsigned long)*8;
    
    
    for(int i=0; i<pRows; i++){
        randomBits[i] = rand()%2; // replace with stronger PRNG.
    }
    
    s = 0;
    enable = 0;
    for(int row = 0; row<pRows; row++){
        if(enable==1)
            break;
        d = 2*d + randomBits[row];
        // Denotes the number of nodes to the right of the node to which we go next.
        for(int col = beginPos[row]; col < pCols; col++) {
            d = d - itos(gel(gel(P, col+1), row+1));
            // We are sort of checking if the hamming distance of all nodes on right of present node. If hamming distance is greater than d, then the current node is a leaf node.
            if(d==-1){
                hit = 1;
                s = col;
                enable = 1;
                break;
            }
        }
        
    }
    //cout<<"Answer sampled : "<<(s*flag)+itos(center)<<endl;
    return (s*flag)+itos(center);
}
