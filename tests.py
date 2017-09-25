import unittest
import Aono

class correctness_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(correctness_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct = Aono.ciphertext(5, self.keys.pk)
        pt = Aono.pari_GEN(5)
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk)
        pt = Aono.pari_GEN([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertEqual(ct.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)

class addition_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(addition_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct_1 = Aono.ciphertext(10, self.keys.pk)
        ct_2 = Aono.ciphertext(5, self.keys.pk)
        pt = Aono.pari_GEN(15)
        ct = ct_1 + ct_2
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct_1 = Aono.ciphertext([1, 2, 3, 4, 5], self.keys.pk)
        ct_2 = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk)
        pt = Aono.pari_GEN([2, 4, 6, 8, 10, 6, 7, 8, 9, 10])
        ct = ct_1 + ct_2
        self.assertEqual(ct.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)

class subtraction_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(subtraction_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct_1 = Aono.ciphertext(10, self.keys.pk)
        ct_2 = Aono.ciphertext(5, self.keys.pk)
        pt = Aono.pari_GEN(5)
        ct = ct_1 - ct_2
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct_1 = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk)
        ct_2 = Aono.ciphertext([1, 2, 3, 4, 5], self.keys.pk)
        pt = Aono.pari_GEN([0, 0, 0, 0, 0, 6, 7, 8, 9, 10])
        ct = ct_1 - ct_2
        self.assertEqual(ct.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)


class miscellaneous(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(miscellaneous, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_ciphertext_array(self):
        ciphertext_array = []
        for i in range(1, 10):
            ciphertext_array.append(Aono.ciphertext(i, self.keys.pk))
        pt = Aono.pari_GEN(5)
        self.assertEqual(ciphertext_array[4].decrypt(self.keys.sk)[0][0], pt)
    
    def test_nested_operations(self):
        ct_1 = Aono.ciphertext(3, self.keys.pk)
        ct_2 = Aono.ciphertext(4, self.keys.pk)
        pt = Aono.pari_GEN(38)
        ct = 2 * ((ct_1 * ct_2) + (ct_1 + ct_2))
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)

class multiplication_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(multiplication_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct_1 = Aono.ciphertext(10, self.keys.pk)
        ct_2 = Aono.ciphertext(5, self.keys.pk)
        pt = Aono.pari_GEN(50)
        ct = ct_1 * ct_2
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct_1 = Aono.ciphertext([1, 2, 3, 4, 5], self.keys.pk)
        ct_2 = Aono.ciphertext([1, 2, 3, 4, 5], self.keys.pk)
        pt = Aono.pari_GEN([[1,2,3,4,5],[2,4,6,8,10],[3,6,9,12,15],[4,8,12,16,20],[5,10,15,20,25]])
        ct = ct_1 * ct_2
        self.assertEqual(ct.decrypt(self.keys.sk).sub_mat_array(0, 5, 0, 5), pt)
    
    def test_single_data_plaintext(self):
        ct = Aono.ciphertext(10, self.keys.pk)
        pt = Aono.pari_GEN(40)
        ct_1 = ct * 4
        ct_2 = 4 * ct
        self.assertEqual(ct_1.decrypt(self.keys.sk)[0][0], pt)
        self.assertEqual(ct_1.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data_plaintext(self):
        ct = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk)
        pt = Aono.pari_GEN([4, 8, 12, 16, 20, 24, 28, 32, 36, 40])
        ct_1 = ct * 4
        ct_2 = 4 * ct
        self.assertEqual(ct_1.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)
        self.assertEqual(ct_2.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)

class keyswitching_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(keyswitching_test, self).__init__(*args, **kwargs)
        self.keysgen = Aono.key_gen()
        self.keys = self.keysgen.generate_key(100, 64, 20, 8, 4, 7)
        self.keys2 = self.keysgen.generate_key(100, 64, 10, 16, 4, 7)
        self.keys3 = self.keysgen.generate_key(100, 64, 25, 8, 4, 7)
        self.updatingkey = Aono.updation_key_gen()
        self.ukey1 = self.updatingkey.generate_key(self.keys, self.keys2)
        self.ukey2 = self.updatingkey.generate_key(self.keys, self.keys3)

    def test_single_data_decrease_security(self):
        ct = Aono.ciphertext(5, self.keys.pk)
        ctupdated = self.ukey1.cipher_switch(ct)
        pt = Aono.pari_GEN(5)
        self.assertEqual(ctupdated.decrypt(self.keys2.sk)[0][0], pt)

    def test_single_data_increase_security(self):
        ct = Aono.ciphertext(5, self.keys.pk)
        ctupdated = self.ukey2.cipher_switch(ct)
        pt = Aono.pari_GEN(5)
        self.assertEqual(ctupdated.decrypt(self.keys3.sk)[0][0], pt)

    def test_multiple_data_increase_security(self):
        ct = Aono.ciphertext([5,4,3,2,1], self.keys.pk)
        ctupdated = self.ukey2.cipher_switch(ct)
        pt = Aono.pari_GEN([5,4,3,2,1])
        self.assertEqual(ctupdated.decrypt(self.keys3.sk).sub_mat_array(0, 5), pt)

    def test_single_data_multiplication(self):
        ct1 = Aono.ciphertext(5, self.keys.pk)
        ct2 = Aono.ciphertext(4, self.keys3.pk)
        ctupdated = self.ukey2.cipher_switch(ct1)
        pt = Aono.pari_GEN(20)
        ct = ctupdated * ct2
        self.assertEqual(ct.decrypt(self.keys3.sk)[0][0], pt)
        ct = 5 * ctupdated
        pt = Aono.pari_GEN(25)
        self.assertEqual(ct.decrypt(self.keys3.sk)[0][0], pt)

    def test_single_data_addition(self):
        ct1 = Aono.ciphertext(5, self.keys.pk)
        ct2 = Aono.ciphertext(4, self.keys3.pk)
        ct3 = Aono.ciphertext(6, self.keys3.pk)
        ctupdated = self.ukey2.cipher_switch(ct1)
        pt = Aono.pari_GEN(9)
        ct = ctupdated + ct2
        self.assertEqual(ct.decrypt(self.keys3.sk)[0][0], pt)
        ctmul = ct2 * ct3
        ct = ctmul + ctupdated
        pt = Aono.pari_GEN(29)
        self.assertEqual(ct.decrypt(self.keys3.sk)[0][0], pt)
        ct = (2 * ctupdated) + ct2
        pt = Aono.pari_GEN(14)
        self.assertEqual(ct.decrypt(self.keys3.sk)[0][0], pt)

    def test_multiple_data_addition(self):
        ct1 = Aono.ciphertext([5,4,3,2], self.keys.pk)
        ct2 = Aono.ciphertext([4,3,2,1], self.keys3.pk)
        ct3 = Aono.ciphertext([6,5,4,3], self.keys3.pk)
        ctupdated = self.ukey2.cipher_switch(ct1)
        pt = Aono.pari_GEN([9,7,5,3])
        ct = ctupdated + ct2
        self.assertEqual(ct.decrypt(self.keys3.sk).sub_mat_array(0, 4), pt)
        ct =  2 * ( 2 * ctupdated + ct2 )
        pt = Aono.pari_GEN([28,22,16,10])
        self.assertEqual(ct.decrypt(self.keys3.sk).sub_mat_array(0, 4), pt)


if __name__ == '__main__':
    unittest.main()

