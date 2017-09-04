import Aono
from Aono import *

pari_init(2000000000, 2)
keys = key_gen()
keys = keys.generate_key(100, 64, 100, 8, 4, 7)

ct = ciphertext(create_message_matrix(4, 64), keys.pk, keys.pk.params)
ct1 = ciphertext(create_message_matrix(4, 64), keys.pk, keys.pk.params)

ct2 = ct + ct1
print_GEN(ct2.decrypt(keys.sk))

