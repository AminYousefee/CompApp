from math import log, exp


class Component:
    required_props_name = ["name", "Tc", "Pc", "omega", "MW", "Tnb", "viscosity25", "Antoine coeffs"]

    def __init__(self):
        self.PHI = 1

    def calc_Dam_Dni(self, fluid):
        n = fluid.n
        am = fluid.am
        ai = self.__a
        ans1 = 2 / n
        ans2 = 0
        for j in range(len(fluid.components)):
            component = fluid.components[j]
            composition = fluid.compositions[j]
            aj = component.a
            ans2 += composition * ((ai * aj) ** 0.5) - am
        return ans1 * ans2

    def calc_Dbm_Dni(self, fluid):
        return (self.b - fluid.bm) / fluid.n

    def calc_Dq_Dni(self, fluid, Dam_Dni, Dbm_Dni):
        term1 = Dam_Dni / fluid.bm
        term2 = fluid.am * Dbm_Dni / (fluid.bm ** 2)
        return term1 - term2

    def calc_DI_Dni(self, EoS, fluid, Dbm_Dni):
        ro = fluid.P / (fluid.T * EoS.R * fluid.z)
        term1 = (ro * EoS.omega) / (1 + EoS.omega * ro * fluid.bm)
        term2 = (ro * EoS.eps) / (1 + EoS.eps * ro * fluid.bm)
        return Dbm_Dni * (term1 - term2)

    def calc_Dz_Dni(self, fluid, Dbm_Dni, EoS, Dq_Dni):
        ro = fluid.P / (fluid.T * EoS.R * fluid.z)
        term1 = ro * Dbm_Dni / (1 - ro * fluid.bm) ** 2
        term2 = ro * fluid.bm * Dq_Dni / ((1 + EoS.eps * ro * fluid.bm) * (1 + EoS.omega * ro * fluid.bm))
        term3 = ro * Dbm_Dni / ((1 + EoS.eps * ro * fluid.bm) * (1 + EoS.omega * ro * fluid.bm))
        term4 = (ro ** 2) * fluid.bm * (EoS.eps ** 2) / (
                ((1 + EoS.eps * ro * fluid.bm) ** 2) * (1 + EoS.omega * ro * fluid.bm))
        term5 = (ro ** 2) * fluid.bm * EoS.omega * Dbm_Dni / (
                (1 + EoS.eps * ro * fluid.bm) * ((1 + EoS.omega * ro * fluid.bm) ** 2))
        term6 = fluid.q * (term3 - term4 - term5)
        term7 = term1 - term2 - term6
        return term7

    def calcu_DGr_Dni(self, fluid, Dz_Dni, Dbm_Dni, EoS, Dq_Dni, DI_Dni):
        ro = fluid.P / (fluid.T * EoS.R * fluid.z)
        term1 = ro * Dbm_Dni * fluid.z / (1 - ro * fluid.bm)
        term2 = 1 / (EoS.omega - EoS.eps)
        I = term2 * log((1 + EoS.omega * ro * fluid.bm) / (1 + EoS.eps * ro * fluid.bm))
        q = EoS.calc_q(fluid.am, fluid.bm, fluid.T)
        term3 = Dz_Dni * log(1 - ro * fluid.bm - I * Dq_Dni - q * DI_Dni)
        ans = EoS.R * fluid.T * (Dz_Dni + term1 - term3)
        return ans

    def calc_DGr_Dni(self, EoS, composition, fluid):
        Dam_Dni = self.calc_Dam_Dni(fluid)
        Dbm_Dni = self.calc_Dbm_Dni(fluid)
        Dq_Dni = self.calc_Dq_Dni(fluid, Dam_Dni, Dbm_Dni)
        DI_Dni = self.calc_DI_Dni(EoS, fluid, Dbm_Dni)
        Dz_Dni = self.calc_Dz_Dni(fluid, Dbm_Dni, EoS, Dq_Dni)
        DGr_Dni = self.calcu_DGr_Dni(fluid, Dz_Dni, Dbm_Dni, EoS, Dq_Dni, DI_Dni)
        return DGr_Dni

    def calc_GbarE_i(self, Gr, ni, DGr_Dni):
        GbarE = Gr + ni * DGr_Dni - self.__Gr
        self.__GbarE = GbarE
        return GbarE

    def calc_gama_i(self, T, EoS):
        self.__gama_i = self.__GbarE / (EoS.R * T)
        return self.__gama_i

    def calc_z_i(self, T, P, EoS):
        a = self.__a
        b = self.__b
        q = EoS.calc_q(a, b, T)
        beta = EoS.calc_beta(b, P, T)
        z = EoS.calc_z_liquid(q, beta)
        self.__z = z
        return z

    def calc_Gr(self, EoS, T, P):
        z = self.__z
        a = self.__a
        b = self.__b
        Gr = EoS.calc_Gr(z, T, P, a, b)
        self.__Gr = Gr
        return Gr

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

    def set_yi_bubble(self, composition, gama, Psat, bubble_P):
        PHI = self.__PHI
        y_i = composition * gama * Psat / (PHI * bubble_P)
        self.__yi_bubble = y_i

    def set_PHI_1(self):
        self.__PHI = 1

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
    def yi_bubble(self):
        return self.__yi_bubble

    @yi_bubble.setter
    def yi_bubble(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def PHI(self):
        return self.__PHI

    @PHI.setter
    def PHI(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def Psat(self):
        return self.__Psat

    @Psat.setter
    def Psat(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def gama_i(self):
        return self.__gama_i

    @gama_i.setter
    def gama_i(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def GbarE(self):
        return self.__GbarE

    @GbarE.setter
    def GbarE(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def Gr(self):
        return self.__Gr

    @Gr.setter
    def Gr(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, value):
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
