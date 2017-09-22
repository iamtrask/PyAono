import Aono
from Aono import *

pari_init(2000000000, 2)
keysgen = key_gen()
keys1 = keysgen.generate_key(100, 64, 40, 8, 4, 7)
message = 11
print "Plaintext message: " + str(message)
ct = ciphertext(create_message_matrix(message, 64), keys1.pk, keys1.pk.params)
ct1 = ciphertext(create_message_matrix(3, 64), keys1.pk, keys1.pk.params)
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

uk = ukgen.generate_key(keys1, keys3)
a = uk.cipher_switch(ct)
b = uk.cipher_switch(ct1)
ct = ciphertext(create_message_matrix(4, 64), keys3.pk, keys3.pk.params)
print "Multiplying rotated-ciphertext with normal-ciphetext"
c = a * ct
print_GEN(get_element(get_element(c.decrypt(keys3.sk), 0), 0))
print "Adding mult-ciphertext with normal-ciphetext"
c1 = c + ct
print_GEN(get_element(get_element(c1.decrypt(keys3.sk), 0), 0))
print "Adding mult-ciphertext with nested normal-ciphetexts"
c2 = (c1 + ct) + ct
print_GEN(get_element(get_element(c2.decrypt(keys3.sk), 0), 0))
print "Adding mult-ciphertext with nested rotated-ciphetexts"
c2 = (c1 + a) + a
print_GEN(get_element(get_element(c2.decrypt(keys3.sk), 0), 0))
print "Adding rotated-ciphertext with normal-ciphetext"
c = a + ct
print_GEN(get_element(get_element(c.decrypt(keys3.sk), 0), 0))
pari_close()




