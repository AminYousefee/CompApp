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

    @property
    def T(self):
        return self.__T

    @T.setter
    def T(self, value):
        raise Exception("you are not allowed to change this property")

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

    def calc_am(self):
        iterator_count = len(self.components)
        sigma = 0
        for j in range(iterator_count):
            for i in range(iterator_count):
                a = ((self.components[i].a * self.components[j].a) ** 0.5)
                sigma += self.compositions[i] * self.compositions[j] * a
        return sigma

    def calc_bm(self):
        sigma = 0
        for i in range(self.components):
            sigma += self.compositions[i] * self.components[i].b
        return sigma

    # every changes execute on this class via main update
    # some properties like Re and Nu is calculated from here
    # fluid is the class which use package property to change and update phases!!!
    # some efficient prop from different phases must be calculate like efficient viscosity
    # this class has a lot of setter and getter and this class is the most concentrated class in this lib
    pass

    def calc_bubble_point(self, T, zi):
        pass

    # calc Pi_sat from....(equation)
    # calc gama_i from...
    #
    #
    # after all calculations which you will write here
    # we have list of composition for each phase and bubble point
    def calc_dew_point(self, T, zi):
        pass
# like bubble point
