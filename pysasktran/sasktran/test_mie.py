import unittest
import sasktran as sk
import numpy as np


class TestMie(unittest.TestCase):
    def test_simple(self):
        mie = sk.MieWiscombe()
        mie.calculate(632.8, 1, 1.5)
