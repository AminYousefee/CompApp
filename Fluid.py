import ProperyContainer


class Fluid(ProperyContainer):
    required_property = dict()

    def __init__(self, T, P, components, compositions):
        if len(components) != len(compositions):
            raise Exception("for each component we need a composition")

    @property
    def Tc(self):
        return self.__Tc

    @Tc.setter
    def Tc(self, value):
        raise Exception("you are not allowed to change this property")

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
