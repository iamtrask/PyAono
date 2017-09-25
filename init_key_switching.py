import Aono
from Aono import *

keysgen = key_gen()
keys1 = keysgen.generate_key(100, 64, 40, 8, 4, 7)
message = 33
print "Plaintext message: " + str(message)
ct = ciphertext(create_message_matrix(message, 64), keys1.pk, keys1.pk.params)

keys2 = keysgen.generate_key(100, 64, 20, 16, 4, 7)
keys3 = keysgen.generate_key(100, 64, 60, 16, 4, 7)
ukgen = updation_key_gen()

print "Reducing security with key switching"
uk = ukgen.generate_key(keys1, keys2)
a = uk.cipher_switch(ct)
print "Decrypted message: "
print_GEN(get_element(get_element(a.decrypt(keys2.sk), 0), 0))

print "Increasing security with key switching"
uk = ukgen.generate_key(keys1, keys3)
a = uk.cipher_switch(ct)
print "Decrypted message: "
print_GEN(get_element(get_element(a.decrypt(keys3.sk), 0), 0))





