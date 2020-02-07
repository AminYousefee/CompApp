import math

from Phase import Phase


class GasPhase(Phase):
    def __init__(self, fluid, compositions, EoS, n):
        super().__init__(fluid, compositions, EoS, n)

    def calc_z(self):
        q = self.EoS.calc_q(self.am, self.bm, self.fluid.T)
        beta = self.EoS.calc_beta(self.bm, self.fluid.P, self.fluid.T)
        z = self.EoS.calc_z_gas(q, beta)




    def calc_kapa_m(self):
        components = self.fluid.components
        kapaj = 0
        kapai = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = self.compositions[i]
                zj = self.compositions[j]
                kapa_i_j = (self.fluid.components[i].kapa * self.fluid.components[j].kapa) ** 0.5
                term1 = zi * zj * (kapa_i_j)
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
                sigma_i = ((self.fluid.components[i].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigma_j = ((self.fluid.components[j].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigma_i_j = (sigma_i * sigma_j) ** 0.5
                term1 = zi * zj * ((sigma_i_j()) ** 3)
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
                sigmai = ((self.fluid.components[i].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigmaj = ((self.fluid.components[j].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                term1 = (zi * zj * mui ** 2 * muj ** 2) / (((sigmai * sigmaj) ** 0.5) ** 3)
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
                sigmai = ((self.fluid.components[i].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigmaj = ((self.fluid.components[j].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                w_i = (self.fluid.components[i].w)
                w_j = (self.fluid.components[j].w)
                w__i_j = (w_i * w_j) / 2
                term1 = (zi * zj * w__i_j * (((sigmai * sigmaj) ** 0.5) ** 3))
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
                sigmai = ((self.fluid.components[i].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigmaj = ((self.fluid.components[j].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                M_ij = (2 * self.fluid.components[i].MW * self.fluid.components[j].MW * 10 ** 3) / (
                        self.fluid.components[i].MW + self.fluid.components[j].MW)
                eps_i = (self.fluid.components[i].Tc) / 1.2593
                eps_j = (self.fluid.components[j].Tc) / 1.2593
                eps_K__i_j = ((eps_i / (1.38 * 10 ** -23)) * (eps_j / (1.38 * 10 ** -23))) ** 0.5
                term1 = (zi * zj * (M_ij ** 0.5) * (((sigmai * sigmaj) ** 0.5) ** 2) * (
                    eps_K__i_j))
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
                eps_i = (self.fluid.components[i].Tc) / 1.2593
                eps_j = (self.fluid.components[j].Tc) / 1.2593
                eps_K__i_j = ((eps_i / (1.38 * 10 ** -23)) * (eps_j / (1.38 * 10 ** -23))) ** 0.5
                sigma_i = ((self.fluid.components[i].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigma_j = ((self.fluid.components[j].Vc * 10 ** 6) ** (1 / 3)) * 0.809
                sigma_i_j = (sigma_i * sigma_j) ** 0.5
                w_i = (self.fluid.components[i].w)
                w_j = (self.fluid.components[j].w)
                w__i_j = (w_i * w_j) / 2

                term1 = (zi * zj * w__i_j * (sigma_i_j ** 3) * eps_K__i_j)
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
        return omega_v
    def calc_y(self):
        ro=self.fluid.ro
        Vc=self.calc_V_c_m()
        y=(ro*Vc)/6
        return y
    def calc_G1(self):
        G1=(1-0.5*self.calc_y())/(1-self.calc_y())**3
        return G1
    def calc_E1(self):
        E1=6.324+50.2412*self.calc_w_m()-51.68*self.calc_mu_4_m()+1189*self.calc_kapa_m()
        return E1
    def calc_E2(self):
        E2=1.21*10**-3-1.154*10**-3*self.calc_w_m()-6.257*10**-3*self.calc_mu_4_m()+0.03728*self.calc_kapa_m()
        return E2
    def calc_E3(self):
        E3=5.283+254.209*self.calc_w_m()-168.48*self.calc_mu_4_m()+3898*self.calc_kapa_m()
        return E3
    def calc_E4(self):
        E4=6.623+38.096*self.calc_w_m()-8.464*self.calc_mu_4_m()+31.42*self.calc_kapa_m()
        return E4
    def calc_E5(self):
        E5=19.745+7.63*self.calc_w_m()-14.354*self.calc_mu_4_m()+31.53*self.calc_kapa_m()
        return E5
    def calc_E6(self):
        E6=-1.9-12.537*self.calc_w_m()+4.985*self.calc_mu_4_m()-18.15*self.calc_kapa_m()
        return E6
    def calc_E7(self):
        E7=24.275+3.45*self.calc_w_m()-11.291*self.calc_mu_4_m()+69.35*self.calc_kapa_m()
        return E7
    def calc_E8(self):
        E8=0.7972+1.117*self.calc_w_m()+0.01235*self.calc_mu_4_m()-4.117*self.calc_kapa_m()
        return E8
    def calc_E9(self):
        E9=-0.2382+0.0677*self.calc_w_m()-0.8163*self.calc_mu_4_m()+4.025*self.calc_kapa_m()
        return E9
    def calc_E10(self):
        E10=0.06863+0.3479*self.calc_w_m()+0.5926*self.calc_mu_4_m()-0.727*self.calc_kapa_m()
        return E10
    def calc_G2(self):
        y=self.calc_y()
        E1=self.calc_E1()
        E2=self.calc_E2()
        E3=self.calc_E3()
        E4=self.calc_E4()
        E5=self.calc_E5()
        G1=self.calc_G1()
        term1=E1*((1-math.e**(-E4*y))/y)
        term2=E2*G1*math.e**(E5*y)
        term3=E3*G1
        term4=E1*E4+E2+E3
        G2=(term1+term2+term3)/term4
        return G2
    def calc_viscosity_star_star(self):
        y=self.calc_y()
        E7=self.calc_E7()
        G2=self.calc_G2()
        E8=self.calc_E8()
        E9=self.calc_E9()
        E10=self.calc_E10()
        T=self.calc_T_star_m()
        term1=E7*y
        term2=E8+E9*T**-1+E10*T**-2
        visco_star_Star=term1*math.e**term2
        return visco_star_Star
    def calc_visco_star(self):
        T=self.calc_T_star_m()
        omega_v=self.calc_omega_v()
        Fc=self.calc_F_c_m()
        G2=self.calc_G2()
        y=self.calc_y()
        E6=self.calc_E6()
        visco_Star_star=self.calc_viscosity_star_star()
        term1=(T**0.5)/omega_v
        term2=(Fc*(G2**-1+E6*y))
        visco_Star=term1*term2+visco_Star_star
        return visco_Star

    def calc_high_pressure_viscosity(self):
        Tc=self.calc_T_c_m()
        M=self.calc_M_m()
        Vc=self.calc_V_c_m()
        visco_star=self.calc_visco_star()
        term1=(M*Tc)**0.5
        term2=Vc**(2/3)
        viscosity=visco_star*((36.344*term1)/term2)
        return viscosity













    def calc_viscosity(self):
        min_high_pressure = 10 * 101325
        viscosity = 0
        if self.fluid.P >= min_high_pressure:
            viscosity = self.calc_high_pressure_viscosity()
        else:
            viscosity = self.calc_low_pressure_viscosity()
        return viscosity
