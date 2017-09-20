import Aono
from Aono import *

pari_init(2000000000, 2)
keysgen = key_gen()
keys1 = keysgen.generate_key(100, 64, 40, 8, 4, 7)


ct = ciphertext(create_message_matrix(26, 64), keys1.pk, keys1.pk.params)

keys2 = keysgen.generate_key(100, 64, 20, 16, 4, 7)

ukgen = updation_key_gen()
uk = ukgen.generate_key(keys1, keys2)
a = uk.cipher_switch(ct)
print_GEN(get_element(get_element(a.decrypt(keys2.sk), 0), 0))
pari_close()




