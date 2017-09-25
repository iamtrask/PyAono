import unittest
import Aono

class correctness_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(correctness_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct = Aono.ciphertext(5, self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN(5)
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertEqual(ct.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)

class addition_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(addition_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct_1 = Aono.ciphertext(10, self.keys.pk, self.keys.pk.params)
        ct_2 = Aono.ciphertext(5, self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN(15)
        ct = ct_1 + ct_2
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct_1 = Aono.ciphertext([1, 2, 3, 4, 5], self.keys.pk, self.keys.pk.params)
        ct_2 = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN([2, 4, 6, 8, 10, 6, 7, 8, 9, 10])
        ct = ct_1 + ct_2
        self.assertEqual(ct.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)

class subtraction_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(subtraction_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct_1 = Aono.ciphertext(10, self.keys.pk, self.keys.pk.params)
        ct_2 = Aono.ciphertext(5, self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN(5)
        ct = ct_1 - ct_2
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data(self):
        ct_1 = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk, self.keys.pk.params)
        ct_2 = Aono.ciphertext([1, 2, 3, 4, 5], self.keys.pk, self.keys.pk.params)
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
            ciphertext_array.append(Aono.ciphertext(i, self.keys.pk, self.keys.pk.params))
        pt = Aono.pari_GEN(5)
        self.assertEqual(ciphertext_array[4].decrypt(self.keys.sk)[0][0], pt)
    
    def test_nested_operations(self):
        ct_1 = Aono.ciphertext(3, self.keys.pk, self.keys.pk.params)
        ct_2 = Aono.ciphertext(4, self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN(38)
        ct = 2 * ((ct_1 * ct_2) + (ct_1 + ct_2))
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)

class multiplication_test(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(multiplication_test, self).__init__(*args, **kwargs)
        keys = Aono.key_gen()
        self.keys = keys.generate_key(100, 64, 40, 8, 4, 7)
    
    def test_single_data(self):
        ct_1 = Aono.ciphertext(10, self.keys.pk, self.keys.pk.params)
        ct_2 = Aono.ciphertext(5, self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN(50)
        ct = ct_1 * ct_2
        self.assertEqual(ct.decrypt(self.keys.sk)[0][0], pt)
    
    def test_single_data_plaintext(self):
        ct = Aono.ciphertext(10, self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN(40)
        ct_1 = ct * 4
        ct_2 = 4 * ct
        self.assertEqual(ct_1.decrypt(self.keys.sk)[0][0], pt)
        self.assertEqual(ct_1.decrypt(self.keys.sk)[0][0], pt)
    
    def test_multiple_data_plaintext(self):
        ct = Aono.ciphertext([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], self.keys.pk, self.keys.pk.params)
        pt = Aono.pari_GEN([4, 8, 12, 16, 20, 24, 28, 32, 36, 40])
        ct_1 = ct * 4
        ct_2 = 4 * ct
        self.assertEqual(ct_1.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)
        self.assertEqual(ct_2.decrypt(self.keys.sk).sub_mat_array(0, 10), pt)

if __name__ == '__main__':
    unittest.main()

