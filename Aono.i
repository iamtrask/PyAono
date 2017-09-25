%module Aono
%{
#include "lib/he/basic.h"
%}
extern void pari_init(size_t parisize, int maxprime);
extern void pari_close();
typedef long *GEN;
%include "lib/he/basic.h"
%include "lib/he/keys.h"
%include "lib/he/utils.h"

%pythoncode%{
    import atexit
    pari_init(2000000000, 2)
    atexit.register(pari_close)
%}

%include "carrays.i"
%array_class(int, intArray);
%extend pari_GEN{
    pari_GEN(PyObject *int_list){
        pari_GEN* result = new pari_GEN();
        if(!PyList_Check(PyList_GetItem( int_list, 0))){
            int *array = NULL;
            int nInts;
            if (PyList_Check( int_list ))
            {
                nInts = PyList_Size( int_list );
                array = (int*) malloc( nInts * sizeof(int) );
                for ( int ii = 0; ii < nInts; ii++ ){
                    PyObject *oo = PyList_GetItem( int_list, ii);
                    if ( PyInt_Check( oo ) )
                        array[ ii ] = ( int ) PyInt_AsLong( oo );
                }
            }
            GEN x;
            x = cgetg(nInts + 1, t_VEC);
            for(int i = 0; i < nInts; i++)
                gel(x, i + 1) = stoi(array[i]);
            result->initialize(x);
            return result;
        }
        else{
            int *array = NULL;
            int **arrofarray = NULL;
            int nInts;
            int nOutInts;
            if (PyList_Check( int_list ))
            {
                nOutInts = PyList_Size( int_list );
                arrofarray = (int**) malloc( nOutInts * sizeof(int*) );
                for ( int jj = 0; jj < nOutInts; jj++ ){
                    nInts = PyList_Size(PyList_GetItem( int_list, jj) );
                    array = (int*) malloc( nInts * sizeof(int) );
                    for ( int ii = 0; ii < nInts; ii++ ){
                        PyObject *oo = PyList_GetItem(PyList_GetItem( int_list, jj), ii);
                        if ( PyInt_Check( oo ) ){
                            array[ ii ] = ( int ) PyInt_AsLong( oo );
                        }
            
                    }
                    arrofarray[jj] = array;
                }
            }
            GEN x;
            x = zeromatcopy(nOutInts, nInts);
            for(int j = 0; j < nOutInts; j++){
                for(int i = 0; i < nInts; i++)
                gel(gel(x, i + 1), j + 1) = stoi(arrofarray[j][i]);
            }
            result->initialize(x);
            return result;
        }
        
    }
    
    char* __str__(){
        return GENtostr(self->value);
    }
    
    pari_GEN __getitem__(int key){
        pari_GEN result;
        result.value = gel(self->value, key + 1);
        return result;
    }
    
    pari_GEN sub_mat_array(int key_1, int key_2){
        pari_GEN result;
        result.value = cgetg(key_2 - key_1 + 1, t_VEC);
        int cnt = 0;
        for(int i = key_1; i < key_2; i++)
            gel(result.value, 1+cnt++) = gel(gel(self->value, i + 1), 1);
        return result;
    }
    
    pari_GEN sub_mat_array(int row_1, int row_2, int col_1, int col_2){
        pari_GEN result;
        result.value = zeromatcopy(row_2 - row_1, col_2 - col_1);
        for(int j = col_1; j < col_2; j++){
            for(int i = row_1; i < row_2; i++){
                gel(gel(result.value, 1+(j-col_1)), 1+(i-row_1)) = gel(gel(self->value, j + 1), i + 1);
            }
        }
        return result;
    }
    
    pari_GEN sub_array(int key_1, int key_2){
        pari_GEN result;
        result.value = cgetg(key_2 - key_1 + 1, t_VEC);
        for(int i = key_1; i < key_2; i++)
            gel(result.value, i + 1) = gel(self->value, i + 1);
        return result;
    }
};

%extend ciphertext{
    ciphertext(PyObject *int_list, public_key* pk){
        ciphertext* result = new ciphertext();
        int *array = NULL;
        int nInts;
        if (PyList_Check( int_list ))
        {
            nInts = PyList_Size( int_list );
            array = (int*) malloc( nInts * sizeof(int) );
            for ( int ii = 0; ii < nInts; ii++ ){
                PyObject *oo = PyList_GetItem( int_list, ii);
                if ( PyInt_Check( oo ) )
                    array[ ii ] = ( int ) PyInt_AsLong( oo );
            }
        }
        pari_GEN pt;
        pt.value = cgetg(nInts + 1, t_VEC);
        for(int i = 0; i < nInts; i++)
            gel(pt.value, i + 1) = stoi(array[i]);
        result->packing_method(pt, pk);
        return result;
    }

    
    ciphertext __mul__(const int pt){
        ciphertext result;
        pari_GEN pt_GEN(pt);
        result.value = plaintext_multiplication(self->value, pt_GEN);
        result.params = self->params;
        result.pk = self->pk;
        return result;
    }
    
    ciphertext __rmul__(const int pt){
        ciphertext result;
        pari_GEN pt_GEN(pt);
        result.value = plaintext_multiplication(self->value, pt_GEN);
        result.params = self->params;
        result.pk = self->pk;
        return result;
    }
};
