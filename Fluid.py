import ProperyContainer


class Fluid(ProperyContainer):
    required_property = dict()

    def __init__(self, T, P, components, compositions):
        if len(components) != len(compositions):
            raise Exception("for each component we need a composition")
        self.__T = T
        self.__P = P
        self.__components = components
        for component in components:
            component.fluid = self
        self.__compositions = compositions

    def calc_a(self, EoS):
        self.calc_alpha(EoS)  # for all component calc alpha
        for component in self.components:
            component.calc_a(EoS)
        iterator_count = len(self.components)
        sigma = 0
        # calc am
        for j in range(iterator_count):
            for i in range(iterator_count):
                other = ((self.components[i].a * self.components[j].a) ** 0.5)
                sigma += self.compositions[i] * self.compositions[j] * other
        self.__am = sigma
        return sigma

    def calc_b(self, EoS):
        iterator_count = len(self.components)
        for component in self.components:
            component.calc_b(EoS)
        sigma = 0
        for i in range(iterator_count):
            sigma += self.components[i].b * self.compositions[i]
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

    # every changes execute on this class via main update
    # some properties like Re and Nu is calculated from here
    # fluid is the class which use package property to change and update phases!!!
    # some efficient prop from different phases must be calculate like efficient viscosity
    # this class has a lot of setter and getter and this class is the most concentrated class in this lib
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

    def calc_Gr_i(self, EoS, ro):
        for i in range(len(self.components)):
            component = self.components[i]
            composition = self.compositions[i]
            component.calc_Gr(EoS, composition, ro, self)

    def calc_bubble_point(self, EoS):
        am = self.calc_a(EoS)
        bm = self.calc_b(EoS)
        z_liquid = EoS.calc_z_liquid(self)
        ro = self.P / (z_liquid * EoS.R * self.T)
        self.calc_Gr_i(EoS, ro)

    # calc Pi_sat from....(equation)
    # calc gama_i from...
    #
    #
    # after all calculations which you will write here
    # we have list of composition for each phase and bubble point
    def calc_dew_point(self, EoS):
        pass
# like bubble point
