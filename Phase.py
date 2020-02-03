import math
from math import log

import ProperyContainer
from EoS import EoS


def calc_H_ig_i(cp_coeffs, T, T0):
    A = cp_coeffs[0]
    B = cp_coeffs[1]
    C = cp_coeffs[2]
    D = cp_coeffs[3]
    H_ig_i = A * (T - T0) + (B / 2) * (T ** 2 - T0 ** 2) + (C / 3) * (T ** 3 - T0 ** 3) - (D / 1) * (
            1 / (T) - (1 / (T0)))
    return H_ig_i


class Phase(ProperyContainer):
    def __init__(self, fluid, compositions, EoS):
        self.__fluid = fluid
        self.__compositions = compositions
        self.__EoS = EoS

    @property
    def fluid(self):
        return self.__fluid

    @fluid.setter
    def fluid(self, fluid):
        raise Exception("you are not allowed to change this property")

    @property
    def EoS(self):
        return self.__EoS

    @EoS.setter
    def EoS(self, EoS):
        raise Exception("you are not allowed to change this property")

    @property
    def compositions(self):
        return self.__compositions

    @compositions.setter
    def compositions(self, compositions):
        raise Exception("you are not allowed to change this property")

    def calc_I(self, EoS, z, bm, fluid):
        P = fluid.P
        T = fluid.T
        ro = P / (T * EoS.R * z)
        term1 = 1 / (EoS.omega - EoS.eps)
        term2 = (EoS.omega * bm * ro + 1) / (1 + EoS.eps * ro * bm)
        term3 = log(term2)
        return term1 * term3

    def calc_Dq_DT(self, am, bm, T, EoS, Dam_DT):
        term1 = Dam_DT / (bm * EoS.R * T)
        term2 = am / (bm * EoS.R * (T ** 2))
        return term1 - term2

    def calc_Dam_DT(self, fluid, compositions):
        T = fluid.T
        components = fluid.components
        sigmaj = 0
        sigmai = 0
        count = len(components)
        for j in range(count):
            for i in range(count):
                zi = compositions[i]
                zj = compositions[j]
                Ai = components[i].a
                Aj = components[j].a
                term1 = zi * zj * ((Ai * Aj) ** 0.5)
                alphaj = components[j].alpha
                alphai = components[i].alpha
                ki = components[i].k
                kj = components[j].k
                Tci = components[i].Tc
                Tcj = components[j].Tc
                term2 = (-ki / (2 * ((T * Tci) ** 0.5)))
                term3 = (-kj / (2 * ((T * Tcj) ** 0.5)))
                sigmai += term1 * (alphaj * term2 + alphai * term3)
            sigmaj += sigmai
        return sigmaj

    def calc_Hr(self, fluid, EoS, z):
        am = fluid.calc_am(EoS, self.compositions)
        bm = fluid.calc_bm(EoS, self.compositions)
        Dam_DT = self.calc_Dam_DT(fluid, self.compositions)
        Dq_DT = self.calc_Dq_DT(am, bm, fluid.T, EoS, Dam_DT)
        I = self.calc_I(EoS, z, bm, fluid)
        return EoS.R * fluid.T * (z - 1 + fluid.T * I * Dq_DT)

    def calc_H_ig(self, fluid):
        H_ig = 0
        component_count = len(fluid.components)
        for i in range(component_count):
            component = fluid.components[i]
            composition = self.compositions[i]
            H_ig_i = calc_H_ig_i(component.Cp_coeff, fluid.T, 298)
            H_ig += composition * H_ig_i
        return H_ig

    def calc_H(self, z):
        EoS = self.__EoS
        fluid = self.fluid
        Hr = self.calc_Hr(fluid, EoS, z)
        H_ig = self.calc_H_ig(fluid)
        return H_ig + Hr

    def calc_entropy_residual(self, fluid, EoS, z):
        bm = fluid.calc_bm(EoS, self.compositions)
        am = fluid.calc_am(EoS, self.compositions)
        Dam_DT = self.calc_Dam_DT(fluid, self.compositions)
        term1 = log(z - EoS.calc_beta(bm, fluid.P, fluid.T))
        term2 = EoS.calc_q(am, bm, fluid.T)
        term3 = fluid.T * (self.calc_Dq_DT(am, bm, fluid.T, EoS, Dam_DT))
        I = self.calc_I(EoS, z, bm, fluid)
        S_R = term1 + (term2 + term3) * I
        return S_R

    def calc_CP_integral_for_S(self, cp_coeffs, T, T0):
        A = cp_coeffs[0]
        B = cp_coeffs[1]
        C = cp_coeffs[2]
        D = cp_coeffs[3]
        cp_intergral_for_s = A * log(T / T0, math.e) + B * (T - T0) + (C / 2) * (T ** 2 - T0 ** 2) - (D / 2) * (
                T ** -2 - T0 ** -2)
        return cp_intergral_for_s

    def calc_p_integral_for_s(self, P, P0):
        P_integral_for_s = EoS.R * (log(P / P0), math.e)
        return P_integral_for_s

    def calc_s_ideal_gas(self, cp_coeffs, T, T0, P, P0):
        S_ig = self.calc_CP_integral_for_S(cp_coeffs, T, T0) - self.calc_p_integral_for_s(P, P0)
        return S_ig

    def calc_S_ig_mixture(self, fluid):
        T = fluid.T
        P = fluid.P
        cp_coeff = fluid.components.Cp_coeffs
        mixture_S_ig = 0
        count = len(fluid.components)
        for i in range(count):
            component = fluid.components[i]
            composition = self.compositions[i]
            S_ig_i = self.calc_s_ideal_gas(component.Cp_coeff, T, 298, P, 1)
            mixture_S_ig += composition * S_ig_i
        return mixture_S_ig

    def calc_mixture_entropy(self, z):
        EoS = self.__EoS
        fluid = self.fluid
        S = self.calc_entropy_residual(fluid, EoS, z) + self.calc_S_ig_mixture(fluid)
        return S
