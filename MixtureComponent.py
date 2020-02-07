import Component


class MixtureComponent(Component):
    def __init__(self, SG, TB):
        self.__SG = SG
        self.__TB = TB

    def calc_MW(self):
        MW = 4.5673 * 10 ** -5 * (self.TB) ** 2.1962 * (self.SG) ** -1.0164
        return MW

    def calc_Tc(self):
        Tc = 24.2787 * (self.TB) ** 0.58848 * (self.SG) ** 0.3596
        return Tc

    def calc_Pc(self):
        Pc = 3.12281 * 10 ** 9 * (self.TB) ** -2.3125 * (self.SG) ** 2.3201
        return Pc

    def calc_Vc(self):
        Vc = 7.5214 * 10 ** -3 * (self.TB) ** 0.2896 * (self.SG) ** -0.7666
        return Vc

    # initial props in this kind of components are less than others and they are calculable
