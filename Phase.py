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
        term3 = log(term2,math.e)
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

    def calc_G_mixture(self, fluid, z):
        G_mixture = self.calc_H(z) - fluid.T * self.calc_mixture_entropy(z)
        return G_mixture

    def calc_U_mixture(self, fluid, z):
        U_mixture = self.calc_H(z) - (z * EoS.R * fluid.T)
        return U_mixture

    def calc_A_mixture(self, fluid, z):
        A_mixture = self.calc_U_mixture(fluid, z) - fluid.T * self.calc_mixture_entropy(z)
        return A_mixture

    def calc_A_for_EOS(self,fluid,EoS):
        am=fluid.calc_am(EoS,self.compositions)
        P=fluid.P
        T=fluid.T
        A=(am*P)/((EoS.R**2)*T**2)
        return A
    def calc_B_for_EOS(self,fluid,EoS):
        bm=fluid.calc_bm(EoS,self.compositions)
        P = fluid.P
        T = fluid.T
        B=(bm*P)/(EoS.R*T)
        return B
    def calc_DB_DT(self,fluid,EoS):
        bm=fluid.calc_bm(EoS,self.compositions)
        P = fluid.P
        T = fluid.T
        DB_DT=-(bm*P)/(EoS.R*T**2)
        return DB_DT
    def calc_DA_DT(self,fluid,EoS):
        P=fluid.P
        T=fluid.T
        R=EoS.R
        am=fluid.calc_am(EoS,self.compositions)
        Dam_DT=self.calc_Dam_DT(fluid,self.compositions)
        term1=P*R**2*T**2*Dam_DT
        term2=2*R**2*T*P*am
        term3=R**4*T**4
        DA_DT=(term1-term2)/term3
        return DA_DT
    def calc_D2am_DT2(self, fluid, compositions):
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
                term4=(0.5*ki / (2 * ((T**1.5) * (Tci ** 0.5))))
                term5 = (0.5 * kj / (2 * ((T ** 1.5) * (Tcj ** 0.5))))
                sigmai += term1 * (alphaj * term4 + alphai * term5+2*(term2*term3))
            sigmaj += sigmai
        return sigmaj
    def calc_D2q_DT2(self,fluid,EoS):
        bm=fluid.calc_bm(EoS,self.compositions)
        am=fluid.calc_am(EoS,self.compositions)
        T=fluid.T
        R=EoS.R
        term1=(bm*R*T)
        term2=(bm*R*T**2)
        term3=self.calc_D2am_DT2(fluid,self.compositions)
        term4=bm*R*self.calc_Dam_DT(fluid,self.compositions)
        term5=self.calc_Dam_DT(fluid,self.compositions)
        term6=2*bm*R*T*am
        D2q_DT2=((term3*term1)-term4)/(term1**2)-((term5*term2)-term6)/(term2**2)
        return D2q_DT2
    def calc_DZ_DT(self,fluid,EoS,z):
        A=self.calc_A_for_EOS(fluid,EoS)
        B=self.calc_B_for_EOS(fluid,EoS)
        DB_DT=self.calc_DB_DT(fluid,EoS)
        DA_DT=self.calc_DA_DT(fluid,EoS)
        term1=(A*DB_DT+B*DA_DT-2*B*DB_DT-3*(B**2)*DB_DT)
        term2=(-DB_DT)*z**2
        term3=z*(DA_DT-2*DB_DT-6*B*DB_DT)
        term4=3*z**2-(1-B)*2*z+(A-2*B-3*B**2)
        DZ_DT=(term1+term2-term3)/term4
        return DZ_DT
    def calc_DI_DT(self,EoS,z,fluid):
        omega=EoS.omega
        eps=EoS.eps
        B=self.calc_B_for_EOS(fluid,EoS)
        DZ_DT=self.calc_DZ_DT(fluid,EoS,z)
        DB_DT=self.calc_DB_DT(fluid,EoS)
        term1=(omega-eps)
        term2=(z+eps*B)
        term3=(z+omega*B)
        term4=(DZ_DT+omega*DB_DT)
        term5=(DZ_DT+eps*DB_DT)
        DI_DT=(1/term1)*((term4/term3)-(term5/term2))
        return DI_DT
    def calc_CP_R(self,z):
        EoS = self.__EoS
        fluid = self.fluid
        T=fluid.T
        R=EoS.R
        bm=fluid.calc_bm(EoS,self.compositions)
        am=fluid.calc_am(EoS,self.compositions)
        Dam_DT=self.calc_Dam_DT(fluid,self.compositions)
        I=self.calc_I(EoS,z,bm,fluid)
        Dq_DT=self.calc_Dq_DT(am,bm,T,EoS,Dam_DT)
        DZ_DT=self.calc_DZ_DT(fluid,EoS,z)
        DI_DT=self.calc_DI_DT(EoS,z,fluid)
        D2q_DT2=self.calc_D2q_DT2(fluid,EoS)
        term1=z-1+T*I*Dq_DT
        term2=DZ_DT+I*Dq_DT+T*I*D2q_DT2+T*DI_DT*Dq_DT
        CP_R=R*term1+R*T*term2
        return CP_R
