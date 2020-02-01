from math import log, exp
class Component:
    required_props_name = ["name", "Tc", "Pc", "omega", "MW", "Tnb", "viscosity25", "Antoine coeffs"]

    def __init__(self):
        pass

    def calc_Gr(self, EoS, composition, ro, fluid):
        a = log(1 - ro * self.__b)
        b = EoS.R * fluid.T * (composition - 1 - composition * a) * composition
        q = fluid.am / (fluid.bm * EoS.R * fluid.T)
        c = 1 / (EoS.omega - EoS.eps)
        d = (1 + EoS.omega *) / ()

    def calc_a(self, EoS):
        a = EoS.a_coeff * 0.45724 * (EoS.R ** 2) * (self.Tc ** 2) / self.Pc
        self.__a = self.__alpha * a
        self.__all_props["a"] = self.__a

    def calc_b(self, EoS):
        self.__b = EoS.b_coeff * EoS.R * self.Tc / self.Pc
        self.__all_props["b"] = self.__b

    def calc_k(self, EoS):
        a = EoS.alpha_coeffs[0]
        b = EoS.alpha_coeffs[1]
        c = EoS.alpha_coeffs[2]
        return a + b * self.omega - c * (self.omega ** 2)

    def calc_alpha(self, k, T):
        Tr = T / self.Tc
        a = 1 - (Tr ** 0.5)
        b = k * a
        c = 1 + b
        self.__alpha = c ** 2
        self.__all_props["alpha"] = self.__alpha
        return self.__alpha

    def calc_Psat(self, T):
        A = self.Antoine_coeffs[0]
        B = self.Antoine_coeffs[1]
        C = self.Antoine_coeffs[2]
        Psat = exp(A - B / (T + C))
        self.__Psat = Psat
        return Psat

    @property
    def Antoine_coeffs(self):
        return self.__Antoine_coeffs

    @Antoine_coeffs.setter
    def Antoine_coeffs(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def alpha(self):
        return self.__alpha

    @alpha.setter
    def alpha(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def all_props(self):
        return self.__all_props

    @all_props.setter
    def all_props(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def Tc(self):
        return self.__Tc

    @Tc.setter
    def Tc(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        raise Exception("you are not allowed to change this property")

    @property
    def Pc(self):
        return self.__Pc

    @Pc.setter
    def Pc(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def omega(self):
        return self.__omega

    @omega.setter
    def omega(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def Tnb(self):
        return self.__Tnb

    @Tnb.setter
    def Tnb(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def MW(self):
        return self.__MW

    @MW.setter
    def MW(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def visco25(self):
        return self.__visco25

    @visco25.setter
    def visco25(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def Cp_coeff(self):
        return self.__Cp_coeff

    @Cp_coeff.setter
    def Cp_coeff(self, value):
        raise Exception("you are not allowed to change this property")

