import math

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

    def calc_w_ij(self):
        count = len(self.fluid.components)
        w_ij = []
        for i in range(count):
            for j in range(count):
                w_i = (self.fluid.components[i].w)
                w_j = (self.fluid.components[j].w)
                w__i_j = (w_i * w_j) / 2
                w_ij.append(w__i_j)
        return w_ij

    def calc_eps_K_ij(self, boltzman):
        count = len(self.fluid.components)
        eps_k_ij = []
        for i in range(count):
            for j in range(count):
                eps_i = (self.fluid.components[i].Tc) / 1.2593
                eps_j = (self.fluid.components[j].Tc) / 1.2593
                eps_K__i_j = ((eps_i / boltzman) * (eps_j / boltzman)) ** 0.5
                eps_k_ij.append(eps_K__i_j)
        return eps_k_ij

    def calc_sigma_ij(self):
        count = len(self.fluid.components)
        sigma_ij = []
        for i in range(count):
            for j in range(count):
                sigma_i = ((self.fluid.components[i].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigma_j = ((self.fluid.components[j].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigma_i_j = (sigma_i * sigma_j) ** 0.5
                sigma_ij.append(sigma_i_j)
        return sigma_ij

    def calc_kapa(self):
        count = len(self.fluid.components)
        kapa_ij = []
        for i in range(count):
            for j in range(count):
                kapa_i_j = (self.fluid.components[i].kapa * self.fluid.components[j].kapa) ** 0.5
                kapa_ij.append(kapa_i_j)
        return kapa_ij

    def calc_M_ij(self):
        count = len(self.fluid.components)
        M_ij = []
        for i in range(count):
            for j in range(count):
                Mij = (2 * self.fluid.components[i].MW * self.fluid.components[j].MW * 10 ** 3) / (
                            self.fluid.components[i].MW + self.fluid.components[j].MW)
                M_ij.append(Mij)
        return M_ij

    def calc_kapa_m(self):
        components = self.fluid.components
        kapaj = 0
        kapai = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                term1 = zi * zj * (self.calc_kapa())
                kapai += term1
            kapaj += kapai
        return kapaj

    def calc_sigma_m(self):
        components = self.fluid.components
        sigmaj = 0
        sigmai = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                term1 = zi * zj * ((self.calc_sigma_ij()) ** 3)
                sigmai += term1
            sigmaj += sigmai
        return sigmaj

    def calc_mu_4_m(self):
        components = self.fluid.components
        mu_4_j = 0
        mu_4_i = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                mui = self.fluid.components[i].mu
                muj = self.fluid.components[j].mu
                term1 = (zi * zj * mui ** 2 * muj ** 2) / ((self.calc_sigma_ij()) ** 3)
                mu_4_i += term1
            mu_4_j += mu_4_i
        mu_4_m = self.calc_sigma_m() ** 3 * mu_4_j
        return mu_4_m

    def calc_w_m(self):
        components = self.fluid.components
        w_j = 0
        w_i = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                term1 = (zi * zj * self.calc_w_ij() * ((self.calc_sigma_ij()) ** 3))
                w_i += term1
            w_j += w_i
        w_m = w_j / self.calc_sigma_m() ** 3
        return w_m

    def calc_M_m(self):
        components = self.fluid.components
        M_j = 0
        M_i = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                term1 = (zi * zj * (self.calc_M_ij() ** 0.5) * ((self.calc_sigma_ij()) ** 2) * (
                    self.calc_eps_K_ij(1.38 * 10 ** -23)))
                M_i += term1
            M_j += M_i
        term2 = self.calc_sigma_m() ** 2 * self.calc_eps_k_m()

        M_m = (M_j / term2) ** 2
        return M_m

    def calc_eps_k_m(self):
        components = self.fluid.components
        eps_k_j = 0
        eps_k_i = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                term1 = (zi * zj * self.calc_w_ij() * ((self.calc_sigma_ij()) ** 3) * (
                    self.calc_eps_K_ij(1.38 * 10 ** -23)))
                eps_k_i += term1
            eps_k_j += eps_k_i
        eps_k_m = eps_k_j / self.calc_sigma_m() ** 3
        return eps_k_m

    def calc_T_star_m(self):
        T_STAR_M = self.fluid.T / (self.calc_eps_k_m())
        return T_STAR_M

    def calc_V_c_m(self):
        V_c_m = (self.calc_sigma_m() / 0.809) ** 3
        return V_c_m

    def calc_T_c_m(self):
        T_c_m = 1.2593 * (self.calc_eps_k_m())
        return T_c_m

    def calc_mu_r_m(self):
        term1 = 131.3 * (self.calc_mu_4_m()) ** 0.25
        term2 = (self.calc_V_c_m() * self.calc_T_c_m()) * 0.5
        mu_r_m = term1 / term2
        return mu_r_m

    def calc_F_c_m(self):
        term1 = 1
        term2 = -0.275 * self.calc_w_m()
        term3 = 0.059035 * (self.calc_mu_r_m()) ** 4
        term4 = self.calc_kapa_m()
        Fcm = term1 + term2 + term3 + term4
        return Fcm

    def calc_low_pressure_viscosity(self):
        Fcm = self.calc_F_c_m()
        Mm = self.calc_M_m()
        T = self.fluid.T
        sigma_m = self.calc_sigma_m()
        omeg_v = self.calc_omega_v()
        viscosity = (26.69 * Fcm * (Mm * T) ** 0.5) / ((sigma_m) ** 2 * omeg_v)
        return viscosity

    def calc_omega_v(self):
        A = 1.15145
        B = 0.14874
        C = 0.52487
        D = 0.7732
        E = 2.16178
        F = 2.43787
        T = self.calc_T_star_m()
        omega_v = (A * (T) ** -B) + C * (math.e) ** (-D * T) + E * (math.e) ** (-F * T)

    def calc_viscosity(self):
        min_high_pressure = 10 * 101325
        viscosity = 0
        if self.fluid.P >= min_high_pressure:
            viscosity = self.calc_high_pressure_viscosity()
        else:
            viscosity = self.calc_low_pressure_viscosity()
        return viscosity
