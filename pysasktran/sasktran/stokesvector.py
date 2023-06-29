import copy
import numpy as np
import sasktranif.sasktranif as skif
from sasktran.exceptions import SasktranError


class StokesVector(object):
    """
    Defines a stokes vector with its associated basis.

    Parameters
    ----------
    stokes : numpy array shape (4,)
        The stokes vector [I, Q, U, V]

    basis : numpy array shape (3,3)
        Coordinate basis the stokes vector is defined in.  basis[0, :] is the propagation direction, basis[1, :] is
        the theta direction, and basis[2, :] is the phi direction.  Directions are specified in ECEF coordinates.
        The basis must be constructed such that basis[1, :] cross basis[2, :] is equal to basis[0, :]

    Raises
    ------
    ValueError
        If the basis is not correctly constructed

    Examples
    --------
    >>> from sasktran import StokesVector
    >>> basis = np.eye(3, 3)
    >>> stokes_rad = [1, 0.1, -0.1, 0]
    >>> stokes_v = StokesVector(stokes_rad, basis)
    >>> print(stokes_v)
    Stokes Vector
    I: 1 Q: 0.1 U: -0.1 V: 0
    Propagation: [ 1.  0.  0.]
    Theta: [ 0.  1.  0.]
    Phi: [ 0.  0.  1.]
    """
    def __init__(self, stokes, basis):
        if not np.isclose(basis[0, :].dot(np.cross(basis[1, :], basis[2, :])), 1):
            raise ValueError('The propagation direction should be equal to the cross product of the theta and phi' +
                             'directions')
        self._stokes = copy.copy(stokes)
        self._basis = copy.copy(basis)

    @classmethod
    def from_skif_object(cls, iskstokesvector: skif.ISKStokesVector):
        """
        Constructs the stokes vector from an already existing `skif.ISKStokesVector` object

        Parameters
        ----------
        iskstokesvector: skif.ISKStokesVector

        Returns
        -------
        StokesVector

        Examples
        --------
        >>> import sasktranif.sasktranif as skif
        >>> from sasktran import StokesVector
        >>> basis = skif.ISKBasisDirection()
        >>> basis.Assign([1, 0, 0], [0, 1, 0], [0, 0, 1])
        >>> stokes_rad = [1, 0.1, -0.1, 0]
        >>> skif_stokes = skif.ISKStokesVector()
        >>> iquv = skif.IQUV()
        >>> iquv.I, iquv.Q, iquv.U, iquv.V = stokes_rad
        >>> skif_stokes.Assign(iquv, basis)
        >>> stokes_v = StokesVector.from_skif_object(skif_stokes)
        >>> print(stokes_v)
        Stokes Vector
        I: 1.0 Q: 0.1 U: -0.1 V: 0.0
        Propagation: [ 1.  0.  0.]
        Theta: [ 0.  1.  0.]
        Phi: [ 0.  0.  1.]
        """
        try:
            stokes = iskstokesvector.Stokes()
            basis = iskstokesvector.Basis()
            basis = np.vstack((basis.Propagation(), basis.Theta(), basis.Phi()))
        except skif._sasktranif.functionfail as e:
            raise SasktranError(e)

        return cls(stokes, basis)

    def __repr__(self):
        representation = ('Stokes Vector\n'
                          'I: {} Q: {} U: {} V: {}\n'
                          'Propagation: {}\n'
                          'Theta: {}\n'
                          'Phi: {}').format(self.I, self.Q, self.U, self.V, self.propagation_direction,
                                            self.theta_direction, self.phi_direction)
        return representation

    def to_new_basis(self, new_basis):
        """
        Converts the stokes vector to a new basis. A rotation matrix between the new basis and old basis is constructed
        and applied to the stokes vector.  Note that this process overrides the old basis

        Parameters
        ----------
        new_basis : numpy array shape (3,3)
            The new basis.  See the class constructor documentation for the format

        Examples
        --------
        >>> from sasktran import StokesVector
        >>> basis = np.eye(3, 3)
        >>> stokes_rad = [1, 0.1, -0.1, 0]
        >>> stokes_v = StokesVector(stokes_rad, basis)
        >>> print(stokes_v)
        Stokes Vector
        I: 1 Q: 0.1 U: -0.1 V: 0
        Propagation: [ 1.  0.  0.]
        Theta: [ 0.  1.  0.]
        Phi: [ 0.  0.  1.]
        >>> new_basis = basis[[0,2,1]]
        >>> new_basis[2] *= -1
        >>> stokes_v.to_new_basis(new_basis)
        >>> print(stokes_v)
        Stokes Vector
        I: 1.0 Q: -0.1 U: 0.1 V: 0.0
        Propagation: [ 1.  0.  0.]
        Theta: [ 0.  0.  1.]
        Phi: [-0. -1. -0.]
        """

        if not np.isclose(new_basis[0, :].dot(self._basis[0, :]), 1):
            raise ValueError('new_basis and current basis should have the same propagation direction')
        if not np.isclose(new_basis[0, :].dot(np.cross(new_basis[1, :], new_basis[2, :])), 1):
            raise ValueError('The propagation direction should be equal to the cross product of the theta and phi' +
                             'directions for the new basis')

        cos_eta = self._basis[1, :].dot(new_basis[1, :])
        sin_eta = -1 * self._basis[1, :].dot(new_basis[2, :])

        cos_two_eta = cos_eta**2 - sin_eta**2
        sin_two_eta = 2.0 * cos_eta * sin_eta

        rot = np.array([[1, 0, 0, 0], [0, cos_two_eta, -sin_two_eta, 0],
                        [0, sin_two_eta, cos_two_eta, 0], [0, 0, 0, 1]])

        self._stokes = rot.dot(self._stokes)
        self._basis = copy.copy(new_basis)

    @property
    def propagation_direction(self):
        """
        Returns the propagation direction of the current basis in ECEF coordinates
        """
        return self._basis[0, :]

    @property
    def theta_direction(self):
        """
        Returns the theta direction of the current basis in ECEF coordinates
        """
        return self._basis[1, :]

    @property
    def phi_direction(self):
        """
        Returns the phi direction of the current basis in ECEF coordinates
        """
        return self._basis[2, :]

    @property
    def I(self):
        return self._stokes[0]

    @property
    def Q(self):
        return self._stokes[1]

    @property
    def U(self):
        return self._stokes[2]

    @property
    def V(self):
        return self._stokes[3]
