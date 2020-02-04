from Phase import Phase


class LiquidPhase(Phase):
    def __init__(self, fluid, compositions, EoS, n):
        super().__init__(fluid, compositions, EoS, n)

    def calc_z(self):
        q = self.EoS.calc_q(self.am, self.bm, self.fluid.T)
        beta = self.EoS.calc_beta(self.bm, self.fluid.P, self.fluid.T)
        z = self.EoS.calc_z_liquid(q, beta)

    def calc_viscosity(self):
        pass
