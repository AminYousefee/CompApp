from Phase import Phase


class GasPhase(Phase):
    def __init__(self, fluid, compositions, EoS, n):
        super().__init__(fluid, compositions, EoS, n)

    def calc_z(self):
        q = self.EoS.calc_q(self.am, self.bm, self.fluid.T)
        beta = self.EoS.calc_beta(self.bm, self.fluid.P, self.fluid.T)
        z = self.EoS.calc_z_gas(q, beta)

    def calc_high_pressure_viscosity(self):
        pass

    def calc_viscosity(self):
        min_high_pressure = 10 * 101325
        viscosity = 0
        if self.fluid.P >= min_high_pressure:
            viscosity = self.calc_high_pressure_viscosity()
        else:
            viscosity = self.calc_low_pressure_viscosity()
        return viscosity
