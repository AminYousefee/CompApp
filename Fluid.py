import ProperyContainer


class Fluid(ProperyContainer):
    required_property = dict()

    def __init__(self, n, T, P, components, compositions):
        if len(components) != len(compositions):
            raise Exception("for each component we need a composition")
        self.__n = n
        self.__T = T
        self.__P = P
        self.__components = components
        for component in components:
            component.fluid = self
        self.__compositions = compositions

    def calc_am(self, EoS, compositions):
        self.calc_alpha(EoS)  # for all component calc alpha
        for component in self.components:
            component.calc_a(EoS)
        iterator_count = len(self.components)
        sigma = 0
        # calc am
        for j in range(iterator_count):
            for i in range(iterator_count):
                other = ((self.components[i].a * self.components[j].a) ** 0.5)
                sigma += compositions[i] * compositions[j] * other
        self.__am = sigma
        return sigma

    def calc_bm(self, EoS, compositions):
        iterator_count = len(self.components)
        for component in self.components:
            component.calc_b(EoS)
        sigma = 0
        for i in range(iterator_count):
            sigma += self.components[i].b * compositions[i]
        self.__bm = sigma
        return sigma

    def calc_alpha(self, EoS):
        for component in self.components:
            k = component.calc_k(EoS)
            component.calc_alpha(k, self.T)

    @property
    def T(self):
        return self.__T

    @T.setter
    def T(self, value):
        self.__T = value

    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, value):
        self.__n = value

    @property
    def P(self):
        return self.__P

    @P.setter
    def P(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def components(self):
        return self.__components

    @components.setter
    def components(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def compositions(self):
        return self.__compositions

    @compositions.setter
    def compositions(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def am(self):
        return self.__am

    @am.setter
    def am(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def bm(self):
        return self.__bm

    @bm.setter
    def bm(self, value):
        raise Exception("you are not allowed to change this property")

    def calc_Gr(self, EoS, z):
        T = self.T
        P = self.P
        a = self.am
        b = self.bm
        Gr = EoS.calc_Gr(z, T, P, a, b)
        return Gr

    def calc_Gr_i(self, EoS, z):
        for component in self.components:
            T = self.T
            P = self.P
            component.calc_Gr(EoS, T, P, z)

    def calc_Psat_i(self):
        for component in self.components:
            component.calc_Psat(self.T)

    def calc_z_liquid(self, EoS):
        a = self.am
        b = self.bm
        T = self.T
        P = self.P
        q = EoS.calc_q(a, b, T)
        beta = EoS.calc_beta(b, P, T)
        z_liq = EoS.calc_z_liquid(q, beta)
        self.__z_liq = z_liq

    def calc_z_i_liq(self, EoS):
        for component in self.components:
            T = self.T
            P = self.P
            component.calc_z_i_liq(T, P, EoS)

    def calc_GbarE_i(self, EoS):
        Gr = self.calc_Gr(EoS, self.__z_liq)
        for i in range(len(self.components)):
            component = self.components[i]
            composition = self.compositions[i]
            ni = self.__n * composition
            DGr_Dni = component.calc_DGr_Dni(EoS, self)
            component.calc_GbarE_i(Gr, ni, DGr_Dni)

    def calc_gama_i(self, EoS):
        for component in self.components:
            component.calc_gama_i(self.T, EoS)

    def calc_bubble_P(self):
        bubble_P = 0
        for i in range(len(self.components)):
            component = self.components[i]
            composition = self.compositions[i]
            gama = component.gama_i
            Psat = component.Psat
            PHI = component.PHI
            bubble_P += composition * gama * Psat / PHI
        return bubble_P

    def calc_bubble_y_i(self, bubble_P):
        self.__bubble_yi = []
        for i in range(self.components):
            component = self.components[i]
            composition = self.compositions[i]
            gama = component.gama_i
            Psat = component.Psat
            y_i = component.set_yi_bubble(composition, gama, Psat, bubble_P)
            self.__bubble_yi.append(y_i)

    def calc_ln_phi_i_hat(self, EoS):
        #  (Gr/ RT) + (n/ RT) * DGr_Dni
        # change am, bm
        bubble_yi =
        amG = self.calc_am(EoS, )
        # change z for mixture
        # change Gr
        new_Gr =
        for component in self.components:
            new_DGr_Dni =

    def calc_bubble_point(self, EoS):
        amL = self.calc_am(EoS, self.compositions)
        bmL = self.calc_bm(EoS, self.compositions)
        self.calc_z_liquid(EoS)
        #  calc log(p sat i) for each component
        self.calc_Psat_i()
        #  calc z for each component(liquid phase)
        self.calc_z_i_liq(EoS)
        #  calc Gr i for each component
        self.calc_Gr_i(EoS, self.__z_liq)
        #  calc Gr for fluid
        #  calc rond(Gr) / rond(ni) for each component
        #  calc Gbar i for each component from 3 up items
        self.calc_GbarE_i(EoS)
        #  calc gama i for each component
        self.calc_gama_i(EoS)
        #  calc P with assumption PHI i = 1
        for component in self.components:
            component.set_PHI_1()
        bubble_P = self.calc_bubble_P()
        self.calc_bubble_y_i(bubble_P)
        self.calc_ln_phi_i_hat()
        #  calc y i for each component
        #

    # calc Pi_sat from....(equation)
    # calc gama_i from...
    #
    #
    # after all calculations which you will write here
    # we have list of composition for each phase and bubble point
    def calc_dew_point(self, EoS):
        pass
# like bubble point
