import numpy as np
from copy import copy


class WignerD:
    def __init__(self, m: int, n: int):
        """
        Class to compute Wigner functions based on the reccurrence relations of Mischenko.  We use the notation
        d^l_{mn} where d^l_{m0} is the associated legendre function.

        Parameters
        ----------
        m : int
        n : int
        """
        self._m = m
        self._n = n

        if self._n >= self._m:
            self._zeta = 1
        else:
            if (self._m - self._n) % 2 == 0:
                self._zeta = 1
            else:
                self._zeta = -1

        self._lmin = max(abs(self._m), abs(self._n))
        self._rec_start_factor = self._recurrence_start_factor()

    def _recurrence_start_factor(self):
        """
        Factor that multiplies the recurrence start that is independent of theta
        """
        num = 2*self._lmin
        den1 = abs(self._m - self._n)
        den2 = abs(self._m + self._n)

        factorial = 1
        for i in range(num, 1, -1):
            factorial *= i
            if i <= den1:
                factorial /= i
            if i <= den2:
                factorial /= i

        otherfactor = 2**(-1*self._lmin)

        return self._zeta * otherfactor * np.sqrt(factorial)

    def _recurrence_start(self, theta: np.array):
        """
        Start factor for the recurrence relations

        Parameters
        ----------
        theta : np.array
            Angle in radians
        """
        x = np.cos(theta)
        xfactor = (1-x)**(abs(self._m - self._n) / 2) * (1+x)**(abs(self._m + self._n) / 2)

        return self._rec_start_factor * xfactor

    def d(self, theta: np.array, l: int):
        """
        Computes d^l_{mn}

        Parameters
        ----------
        theta : np.array
            Angles in radians
        l : int

        """
        if l < self._lmin:
            return 0

        x = np.cos(theta)
        val_l = self._recurrence_start(theta)
        val_lm1 = np.zeros_like(theta)

        if self._n == 0:
            for lidx in range(self._lmin + 1, l+1):
                multiplier = 1 / (np.sqrt(lidx*lidx - self._m*self._m) * lidx)
                curfactor = (2*lidx - 1) * (lidx * x)
                priorfactor = lidx * np.sqrt((lidx-1)*(lidx-1) - self._m*self._m)

                temp = copy(val_l)
                val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1)
                val_lm1 = temp
        else:
            for lidx in range(self._lmin + 1, l+1):
                multiplier = 1 / ((lidx-1) * np.sqrt(lidx*lidx - self._m*self._m) * np.sqrt(lidx*lidx - self._n*self._n))
                curfactor = (2*lidx - 1) * (lidx * (lidx-1) * x - self._n * self._m)
                priorfactor = lidx * np.sqrt((lidx-1) * (lidx-1) - self._m*self._m) * np.sqrt((lidx-1) * (lidx-1) - self._n*self._n)

                temp = copy(val_l)
                val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1)
                val_lm1 = temp
        return val_l
