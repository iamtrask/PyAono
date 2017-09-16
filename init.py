import Aono
from Aono import *

pari_init(2000000000, 2)
keys = key_gen()
keys = keys.generate_key(100, 64, 100, 8, 4, 7)

ct1 = ciphertext(create_message_matrix(4, 64), keys.pk, keys.pk.params)
ct2 = ciphertext(create_message_matrix(5, 64), keys.pk, keys.pk.params)
ct3 = ciphertext(create_message_matrix(6, 64), keys.pk, keys.pk.params)

run = True

while run:

    x = input("Press (1, 2, 3, 4, 5 or 6) for next demo: ")

    if x == 1:
        print("\n-----------------------------Addition(4+5)---------------------------\n")
        ctres = (ct1 + ct2)
        print_GEN(get_element(get_element(ctres.decrypt(keys.sk), 0), 0))

    if x == 2:
        print("\n-----------------------------Subtraction(5-4)---------------------------\n")
        ctres = (ct2 - ct1)
        print_GEN(get_element(get_element(ctres.decrypt(keys.sk), 0), 0))

    if x == 3:
        print("\n-----------------------------Multiplication(4*5)---------------------------\n")
        ctres = (ct1 * ct2)
        print_GEN(get_element(get_element(ctres.decrypt(keys.sk), 0), 0))

    if x is 4:
        print("\n------------------Nested Addition and Multiplication(((4+5)+5)*6)---------------------------\n")
        ctres = ((ct1 + ct2) + ct2 ) * ct3
        # first get_element gives the first column and to get the first element, we use another get_element.
        print_GEN(get_element(get_element(ctres.decrypt(keys.sk), 0), 0))

    if x == 5:
        print("\n------------------Example ciphertext---------------------------\n")
        print_GEN(see_ciphertext(ctres.value, 0))

    if x == 6:
        run = False

pari_close()




