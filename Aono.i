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
