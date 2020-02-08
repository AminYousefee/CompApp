from math import *

from Phase import Phase


class LiquidPhase(Phase):
    def __init__(self, fluid, compositions, EoS, n):
        super().__init__(fluid, compositions, EoS, n)

    def calc_z(self):
        q = self.EoS.calc_q(self.am, self.bm, self.fluid.T)
        beta = self.EoS.calc_beta(self.bm, self.fluid.P, self.fluid.T)
        z = self.EoS.calc_z_liquid(q, beta)

    def calc_H1(self, T, Tc):
        Tr = T / Tc
        A = 0.33593
        B = 0.33953
        C = 1.51941
        D = 0.02512
        E = 1011422
        return A - B * Tr + C * (Tr ** 2) - D * (Tr ** 3) + E * (Tr ** 4)

    def calc_H2(self, T, Tc):
        Tr = T / Tc
        return 0.29067 - 0.09045 * Tr - 0.04842 * (Tr ** 2)

    def calc_F(self, T, H1, H2, omega):
        return H1 * (1 - omega * H2)

    def calc_V_T(self, T, TR, VR, component):
        Tc = component.Tc
        H1 = self.calc_H1(T, Tc)
        H2 = self.calc_H2(T, Tc)
        H1R = self.calc_H1(TR, Tc)
        H2R = self.calc_H2(TR, Tc)
        omega = component.omega
        return self.calc_F(T, H1, H2, omega) / (self.calc_F(TR, H1R, H2R, omega)) * VR

    def calc_V0(self, component):
        omega = component.omega
        Tc = component.Tc
        Tfp = component.Tfp
        Vm = component.Vm
        term1 = 0.0085 * omega * Tc - 2.02
        term2 = Vm / (0.342 * (Tfp / Tc) + 0.894)
        return term1 + term2

    def calc_E(self, component, T):
        MW = component.MW
        Tc = component.Tc
        Tfp = component.Tfp
        Pc = component.Pc
        Vc = component.Vc
        return -1.12 * Vc / (1294 + 0.1 * MW - 0.23 * Pc + 0.0424 * Tfp - 11.58 * (Tfp / T))

    def calc_visc_liq(self, component, T, TR, VR):
        V0 = self.calc_V0(component)
        V = self.calc_V_T(T, TR, VR, component)
        E = self.calc_E(component, T)
        return V0 / (E * V - V0)

    def calc_Vcij(self, Vci, Vcj):
        return ((Vci ** (1 / 3) + Vcj ** (1 / 3)) ** 3) / 8

    def calc_Vcm(self):
        sigma = 0
        count = len(self.compositions)
        for i in range(count):
            for j in range(count):
                xi = self.compositions[i]
                xj = self.compositions[j]
                Vci = self.fluid.components[i].Vc
                Vcj = self.fluid.components[j].Vc
                Vcij = self.calc_Vcij(Vci, Vcj)
                sigma += xi * xj * Vcij
        return sigma

    def calc_omega_m(self):
        count = len(self.compositions)
        sigma = 0
        for i in range(count):
            xi = self.compositions[i]
            omega_i = self.fluid.components[i].omega
            sigma += xi * omega_i
        return sigma

    def calc_MW_m(self):
        count = len(self.compositions)
        sigma = 0
        for i in range(count):
            xi = self.compositions[i]
            MW_i = self.fluid.components[i].MW
            sigma += xi * MW_i
        return sigma

    def calc_Tcm(self):
        sigma = 0
        count = len(self.compositions)
        for i in range(count):
            for j in range(count):
                xi = self.compositions[i]
                xj = self.compositions[j]
                Tci = self.fluid.components[i].Tc
                Tcj = self.fluid.components[j].Tc
                Vci = self.fluid.components[i].Vc
                Vcj = self.fluid.components[j].Vc
                factor = (Tci * Tcj * Vci * Vcj) ** 0.5
                sigma += xi * xj * factor
        Vcm = self.calc_Vcm()
        return sigma / Vcm

    def calc_viscosity(self):
        VR = float(input())
        if len(self.compositions) == 1:
            component = self.fluid.components[0]
            return self.calc_visc_liq(component, self.fluid.T, 298, VR)
        component1 = self.fluid.components[0]
        component2 = self.fluid.components[1]
        Tc1 = component1.Tc
        Tc2 = component2.Tc
        Tcm = self.calc_Tcm()
        T1 = self.fluid.T * Tc1 / Tcm
        T2 = self.fluid.T * Tc2 / Tcm
        visc1 = self.calc_visc_liq(component1, T1, 298, VR)
        visc2 = self.calc_visc_liq(component2, T2, 298, VR)
        eps1 = component1.Vc ** (2 / 3) / (Tc1 * component1.MW) ** 0.5
        eps2 = component2.Vc ** (2 / 3) / (Tc2 * component2.MW) ** 0.5
        omega_m = self.calc_omega_m()
        omega1 = component1.omega
        omega2 = component2.omega
        Vcm = self.calc_Vcm()
        MW_m = self.calc_MW_m()
        eps_m = Vcm ** (2 / 3) / (Tcm * MW_m) ** 0.5
        term1 = log(visc1 * eps1)
        term2 = log(visc2 * eps2) - term1
        term3 = (omega_m - omega1) / (omega2 - omega1)
        term4 = term1 + term2 * term3
        return exp(term4) / eps_m
