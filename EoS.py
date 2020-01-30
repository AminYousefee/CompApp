class EoS:
    R = 8.3145

    def __init__(self, a_coeff, b_coeff, alpha_coeffs, eps, omega, ):
        self.__a_coeff = a_coeff
        self.__b_coeff = b_coeff
        self.__alpha_coeffs = alpha_coeffs
        self.__eps = eps
        self.__omega = omega

    @property
    def a_coeff(self):
        return self.__a_coeff

    @a_coeff.setter
    def a_coeff(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def b_coeff(self):
        return self.__b_coeff

    @b_coeff.setter
    def b_coeff(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def alpha_coeffs(self):
        return self.__alpha_coeffs

    @alpha_coeffs.setter
    def alpha_coeffs(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def eps(self):
        return self.__eps

    @eps.setter
    def eps(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def omega(self):
        return self.__omega

    @omega.setter
    def omega(self, value):
        raise Exception("you are not allowed to change this property")

    def a_calculator(self, component):
        return self.a_coeff * (EoS.R ** 2) * (component.Tc ** 2) / component.Pc

    def b_calculator(self, component):
        return self.b_coeff * EoS.R * component.Tc / component.P4

    def k_calculator(self, component):
        a = self.alpha_coeffs[0]
        b = self.alpha_coeffs[1]
        c = self.alpha_coeffs[2]
        return a + b * component.omega - c * (component.omega ** 2)

    def alpha_calculator(self, component):
        T = component.fluid.T
        Tr = T / component.Tc
        a = 1 - Tr ** 0.5
        b = component.k * a
        c = 1 + b
        return c ** 2
    # now we need to calc alpha for each component
    # we need tempreture but temp is belonged to fluid and here we don't access Fluid
    #
