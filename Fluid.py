from math import exp

from GasPhase import GasPhase
from LiquidPhase import LiquidPhase


class Fluid:

    def __init__(self, n, T, P, components, compositions, EoS, container):
        self.__phase_list = []
        if len(components) == 1:
            self.calc_Psat_i(self.T)
            self.calc_alpha(EoS)
            component = components[0]
            Psat = component.Psat
            if P > Psat:
                # liquid
                liquid = LiquidPhase(self, [1], EoS, n)
                self.__phase_list.append(liquid)
            elif P == Psat:
                total_v = container.V
                z_gas = self.calc_z_gas(EoS, component.a, component.b, T, P)
                z_liq = self.calc_z_liquid(EoS, component.a, component.b, T, P)
                v_gas = (z_gas * EoS.R * T) / P
                v_liq = (z_liq * EoS.R * T) / P
                V = (total_v - n * v_liq) / (v_gas - v_liq)
                L = n - V
                liquid = LiquidPhase(self, [1], EoS, L)
                gas = GasPhase(self, [1], EoS, V)
                self.__phase_list.append(gas)
                self.__phase_list.append(liquid)
            else:
                gas = GasPhase(self, [1], EoS, n)
                # gas
        else:
            if len(components) != len(compositions):
                raise Exception("for each component we need a composition")
            self.__n = n
            self.__T = T
            self.__P = P
            self.__components = components
            self.__EoS = EoS
            self.__compositions = compositions
            self.__bubble_dew_pressure_acc = 1
            self.update_all(self.__bubble_dew_pressure_acc)

    def update_all(self, bubble_dew_pressure_acc):
        self.__phase_list = self.flash_calc(self.__EoS, bubble_dew_pressure_acc)

    def update(self, new_condition):
        for prop in new_condition:
            if prop == "T":
                self.__T = new_condition["T"]
            if prop == "P":
                self.__P = new_condition["P"]
            if prop == "composition":
                if len(self.compositions) != len(new_condition["composition"]):
                    raise Exception("number of compositions is not match with component number")
                self.__compositions = new_condition["composition"]
            self.update_all(self.__bubble_dew_pressure_acc)

    def calc_am(self, EoS, compositions):
        self.calc_alpha(EoS)  # for all component calc alpha
        for component in self.components:
            component.calc_a(EoS)
        iterator_count = len(self.components)
        sigma = 0
        # calc am
        for j in range(iterator_count):
            for i in range(iterator_count):
                other = ((self.components[i].a * self.components[j].a) ** 0.5)
                sigma += compositions[i] * compositions[j] * other
        return sigma

    def calc_bm(self, EoS, compositions):
        iterator_count = len(self.components)
        for component in self.components:
            component.calc_b(EoS)
        sigma = 0
        for i in range(iterator_count):
            sigma += self.components[i].b * compositions[i]
        return sigma

    def calc_alpha(self, EoS):
        for component in self.components:
            k = component.calc_k(EoS)
            component.calc_alpha(k, self.T)

    @property
    def T(self):
        return self.__T

    @T.setter
    def T(self, value):
        self.__T = value

    @property
    def phase_list(self):
        return self.__phase_list

    @phase_list.setter
    def phase_list(self, value):
        self.__phase_list = value

    @property
    def EoS(self):
        return self.__EoS

    @EoS.setter
    def EoS(self, EoS):
        self.__EoS = EoS

    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, value):
        self.__n = value

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

    @property
    def am(self):
        return self.__am

    @am.setter
    def am(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def bm(self):
        return self.__bm

    @bm.setter
    def bm(self, value):
        raise Exception("you are not allowed to change this property")

    @property
    def phase_list(self):
        return self.__phase_list

    @phase_list.setter
    def phase_list(self, value):
        raise Exception("you are not allowed to change this property")

    def calc_Gr(self, EoS, z):
        T = self.T
        P = self.P
        a = self.am
        b = self.bm
        Gr = EoS.calc_Gr(z, T, P, a, b)
        return Gr

    def calc_Gr_i(self, EoS, list_zi, T, P):
        list_Gr_i = []
        for i in range(len(self.components)):
            component = self.components[i]
            zi = list_zi[i]
            Gr_i = component.calc_Gr(EoS, T, P, zi)
            list_Gr_i.append(Gr_i)
        return list_Gr_i

    def calc_Psat_i(self, T):
        for component in self.components:
            component.calc_Psat(T)

    def calc_z_gas(self, EoS, a, b, T, P):
        q = EoS.calc_q(a, b, T)
        beta = EoS.calc_beta(b, P, T)
        z_gas = EoS.calc_z_gas(q, beta)
        return z_gas

    def calc_z_liquid(self, EoS, a, b, T, P):
        q = EoS.calc_q(a, b, T)
        beta = EoS.calc_beta(b, P, T)
        z_liq = EoS.calc_z_liquid(q, beta, b)
        return z_liq

    def calc_z_i_gas(self, EoS, T, P):
        list_z_i_gas = []
        for component in self.components:
            a = component.a
            b = component.b
            q = EoS.calc_q(a, b, T)
            beta = EoS.calc_beta(b, P, T)
            zi_gas = EoS.calc_z_liquid(q, beta, b)
            list_z_i_gas.append(zi_gas)
        return list_z_i_gas

    def calc_z_i_liq(self, EoS, T, P):
        list_z_i_liq = []
        for component in self.components:
            a = component.a
            b = component.b
            q = EoS.calc_q(a, b, T)
            beta = EoS.calc_beta(b, P, T)
            zi_liq = EoS.calc_z_liquid(q, beta, b)
            list_z_i_liq.append(zi_liq)
        return list_z_i_liq

    def calc_GbarE_i(self, Gr, list_Gr_i, list_DGri_Dni):
        n = self.__n
        list_GbarE_i = []
        for i in range(len(self.components)):
            DGri_Dni = list_DGri_Dni[i]
            Gr_i = list_Gr_i[i]
            GbarE_i = Gr + n * DGri_Dni - Gr_i
            list_GbarE_i.append(GbarE_i)
        return list_GbarE_i

    def calc_gama_i(self, EoS, T, list_GbarE_i):
        R = EoS.R
        list_gama_i = []
        for i in range(len(self.components)):
            GbarE_i = list_GbarE_i[i]
            ln_gama_i = GbarE_i / (R * T)
            gama_i = exp(ln_gama_i)
            list_gama_i.append(gama_i)
        return list_gama_i

    def calc_bubble_P(self, composition, list_gama_i, PHI_i_list):
        bubble_P = 0
        for i in range(len(self.components)):
            component = self.components[i]
            gama = list_gama_i[i]
            Psat = component.Psat
            PHI = PHI_i_list[i]
            bubble_P += composition * gama * Psat / PHI
        return bubble_P

    def calc_bubble_y_i(self, bubble_P, composition, list_gama_i, PHI_i_list):
        list_bubble_yi = []
        for i in range(self.components):
            component = self.components[i]
            gama = list_gama_i[i]
            Psat = component.Psat
            PHI = PHI_i_list[i]
            y_i = component.set_yi_bubble(composition, gama, Psat, bubble_P, PHI)
            list_bubble_yi.append(y_i)
        return list_bubble_yi

    def calc_Dam_Dni(self, am, composition):
        list_Dam_Dni = []
        for component in self.components:
            Dam_Dni = component.calc_Dam_Dni(self, am, composition)
            list_Dam_Dni.append(Dam_Dni)
        return list_Dam_Dni

    def calc_Dbm_Dni(self, bm):
        list_Dbm_Dni = []
        for component in self.components:
            Dbm_Dni = component.calc_Dbm_Dni(self, bm)
            list_Dbm_Dni.append(Dbm_Dni)
        return list_Dbm_Dni

    def calc_Dq_Dni(self, list_Dam_Dni, list_Dbm_Dni, am, bm):
        Dq_list = []
        for i in range(len(self.components)):
            component = self.components[i]
            Dam_Dni = list_Dam_Dni[i]
            Dbm_Dni = list_Dbm_Dni[i]
            Dq_Dni = component.calc_Dq_Dni(Dam_Dni, Dbm_Dni, am, bm)
            Dq_list.append(Dq_Dni)
        return Dq_list

    def calc_DI_Dni(self, EoS, T, z, list_Dbm_Dni, bm, P):
        list_DI = []
        for i in range(len(self.components)):
            component = self.components[i]
            Dbm_Dni = list_Dbm_Dni[i]
            DI_Dni = component.calc_DI_Dni(EoS, T, P, z, bm, Dbm_Dni)
            list_DI.append(DI_Dni)
        return list_DI

    def calc_Dz_Dni(self, list_Dbm_Dni, EoS, list_Dq_Dni, a, b, T, P, z):
        list_Dz = []
        for i in range(self.components):
            component = self.components[i]
            Dbm_Dni = list_Dbm_Dni[i]
            Dq_Dni = list_Dq_Dni[i]
            Dz = component.calc_Dz_Dni(Dbm_Dni, EoS, Dq_Dni, a, b, T, P, z)
            list_Dz.append(Dz)
        return list_Dz

    def calcu_DGri_Dni(self, a, b, z, T, P, composition, EoS):
        list_Dam_Dni = self.calc_Dam_Dni(self, composition)
        list_Dbm_Dni = self.calc_Dbm_Dni(b)
        list_Dq_Dni = self.calc_Dq_Dni(list_Dam_Dni, list_Dbm_Dni, a, b)
        list_DI_Dni = self.calc_DI_Dni(EoS, T, z, list_Dbm_Dni, b, P)
        list_Dz_Dni = self.calc_Dz_Dni(list_Dbm_Dni, EoS, list_Dq_Dni, a, b, T, P, z)
        list_DGri_Dni = []
        for i in range(len(self.components)):
            Dbm = list_Dbm_Dni[i]
            Dq = list_Dq_Dni[i]
            DI = list_DI_Dni[i]
            Dz = list_Dz_Dni[i]
            component = self.components[i]
            DGri_Dni = component.calcu_DGr_Dni(Dz, Dbm, EoS, Dq, DI, T, P, z, a, b)
            list_DGri_Dni.append(DGri_Dni)
        return list_DGri_Dni

    def calc_DGri_Dni(self, a, b, z, composition, EoS):
        list_DGri_Dni = []
        for i in range(len(self.components)):
            DGri_Dni = self.calcu_DGri_Dni(a, b, z, self.T, self.P, composition, EoS)
            list_DGri_Dni.append(DGri_Dni)
        return list_DGri_Dni

    def calc_phi_i_hat(self, Gr_gas, T, EoS, n, list_DGr_Dni_gas):
        list_phi_i = []
        for i in range(len(self.components)):
            DGr_Dni = list_DGr_Dni_gas[i]
            ln_phi = Gr_gas / (EoS.R * T) + n * DGr_Dni / (EoS.R * T)
            phi = exp(ln_phi)
            list_phi_i.append(phi)
        return list_phi_i

    def calc_ln_phi_i_hat_sat(self, Gr_gas, T, EoS, n, am, bm, z_gas):
        list_DGri_Dni = []
        for i in range(len(self.components)):
            yi_list = [0] * len(self.components)
            yi_list[i] = 1
            list_DGri_Dni_gas2 = self.calc_DGri_Dni(am, bm, z_gas, yi_list, EoS)
            DGri_Dni_gas = list_DGri_Dni_gas2[i]
            list_DGri_Dni.append(DGri_Dni_gas)
        list_phi_i_hat_sat = self.calc_phi_i_hat(Gr_gas, T, EoS, n, list_DGri_Dni)
        return list_phi_i_hat_sat

    def calc_PHI(self, list_phi_i_hat, list_phi_i_hat_sat):
        list_PHI = []
        for i in range(len(self.components)):
            PHI = list_phi_i_hat[i] / list_phi_i_hat_sat
            list_PHI.append(PHI)
        return list_PHI

    def calc_bubble_point(self, EoS, acc):
        amL = self.calc_am(EoS, self.compositions)
        bmL = self.calc_bm(EoS, self.compositions)
        z_liq = self.calc_z_liquid(EoS, amL, bmL, self.T, self.P)
        self.calc_Psat_i(self.T)
        list_zi_liq = self.calc_z_i_liq(EoS, self.T, self.P)
        Gr_liq = self.calc_Gr(EoS, z_liq)
        list_Gr_i_liq = self.calc_Gr_i(EoS, list_zi_liq, self.T, self.P)
        list_DGri_Dni_liq = self.calc_DGri_Dni(amL, bmL, z_liq, self.compositions, EoS)
        list_GbarE_i = self.calc_GbarE_i(Gr_liq, list_Gr_i_liq, list_DGri_Dni_liq)

        list_gama_i = self.calc_gama_i(EoS, self.T, list_GbarE_i)

        PHI_list = [1] * len(self.components)
        #####################
        #####################
        # old(initial)
        bubble_P_old = self.calc_bubble_P(self.compositions, list_gama_i, PHI_list)
        list_bubble_yi = self.calc_bubble_y_i(bubble_P_old, self.compositions, list_gama_i, PHI_list)
        amV = self.calc_am(EoS, list_bubble_yi)
        bmV = self.calc_bm(EoS, list_bubble_yi)
        z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, bubble_P_old)
        Gr_gas = self.calc_Gr(EoS, z_gas)
        # list_zi_gas = self.calc_z_i_gas(EoS, self.T, bubble_P)
        # list_Gr_i_gas = self.calc_Gr_i(EoS, list_zi_gas, self.T, bubble_P)
        list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, list_bubble_yi, EoS)
        list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, self.n, list_DGri_Dni_gas1)
        list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, self.n, amV, bmV, z_gas)
        list_PHI_i = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)
        ##########################
        ##########################
        # new
        bubble_P_new = self.calc_bubble_P(self.compositions, list_gama_i, list_PHI_i)
        list_bubble_yi_new = self.calc_bubble_y_i(bubble_P_new, self.compositions, list_gama_i, list_PHI_i)
        amV = self.calc_am(EoS, list_bubble_yi_new)
        bmV = self.calc_bm(EoS, list_bubble_yi_new)
        z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, bubble_P_new)
        Gr_gas = self.calc_Gr(EoS, z_gas)
        list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, list_bubble_yi_new, EoS)
        list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, self.n, list_DGri_Dni_gas1)
        list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, self.n, amV, bmV, z_gas)
        list_PHI_i = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)
        ######################
        ######################
        while abs(bubble_P_new - bubble_P_old) > acc:
            bubble_P_old = bubble_P_new
            bubble_P_new = self.calc_bubble_P(self.compositions, list_gama_i, list_PHI_i)
            list_bubble_yi_new = self.calc_bubble_y_i(bubble_P_new, self.compositions, list_gama_i, list_PHI_i)
            amV = self.calc_am(EoS, list_bubble_yi_new)
            bmV = self.calc_bm(EoS, list_bubble_yi_new)
            z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, bubble_P_new)
            Gr_gas = self.calc_Gr(EoS, z_gas)
            list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, list_bubble_yi_new, EoS)
            list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, self.n, list_DGri_Dni_gas1)
            list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, self.n, amV, bmV, z_gas)
            list_PHI_i = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)
            ########
            ########
        return [bubble_P_new, list_PHI_i, list_gama_i, list_bubble_yi_new]

    def calc_Dew_P(self, composition, PHI_list, list_gama):
        sigma = 0
        for i in range(len(self.components)):
            component = self.components[i]
            Psat = component.Psat
            yi = composition[i]
            PHI_i = PHI_list[i]
            gama_i = list_gama[i]
            sigma += ((yi * PHI_i) / (gama_i * Psat))
        return 1 / sigma

    def calc_x(self, yi_list, PHI_list, dew_P, gama_list):
        list_xi = []
        for i in range(len(self.components)):
            yi = yi_list[i]
            PHI = PHI_list[i]
            gama = gama_list[i]
            component = self.components[i]
            Psat = component.Psat
            xi = (yi * PHI * dew_P) / (gama * Psat)
            list_xi.append(xi)
        return list_xi

    def calc_gama(self, EoS, list_GbarE_i):
        T = self.T
        R = EoS.R
        list_gama = []
        for i in range(len(self.components)):
            component = self.components[i]
            GbarE_i = list_GbarE_i[i]
            gama = GbarE_i / (R * T)
            list_gama.append(gama)
        return list_gama

    def calc_dew_point(self, EoS, bubble_dew_pressure_acc):
        count = len(self.components)
        yi_list = self.compositions
        ############
        ############
        # initial(old)
        list_PHI_old = [1] * count
        list_gama_old = [1] * count

        Dew_P_old = self.calc_Dew_P(yi_list, list_PHI_old, list_gama_old)
        list_x_old = self.calc_x(yi_list, list_PHI_old, Dew_P_old, list_gama_old)

        amL = self.calc_am(EoS, list_x_old)
        bmL = self.calc_bm(EoS, list_x_old)
        z_liq = self.calc_z_liquid(EoS, amL, bmL, self.T, Dew_P_old)
        list_zi_liq = self.calc_z_i_liq(EoS, self.T, Dew_P_old)
        Gr_liq = self.calc_Gr(EoS, z_liq)
        list_Gr_i_liq = self.calc_Gr_i(EoS, list_zi_liq, self.T, Dew_P_old)
        list_DGri_Dni_liq = self.calc_DGri_Dni(amL, bmL, z_liq, list_x_old, EoS)
        list_GbarE_i = self.calc_GbarE_i(Gr_liq, list_Gr_i_liq, list_DGri_Dni_liq)
        list_gama = self.calc_gama(EoS, list_GbarE_i)

        Dew_P_old = self.calc_Dew_P(yi_list, list_PHI_old, list_gama)
        # list_x_old = self.calc_x(yi_list, list_PHI_old, Dew_P_old, list_gama)

        amV = self.calc_am(EoS, yi_list)
        bmV = self.calc_bm(EoS, yi_list)
        z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, self.P)
        Gr_gas = self.calc_Gr(EoS, z_gas)
        # list_zi_gas = self.calc_z_i_gas(EoS, self.T, bubble_P)
        # list_Gr_i_gas = self.calc_Gr_i(EoS, list_zi_gas, self.T, bubble_P)
        list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, yi_list, EoS)
        list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, self.n, list_DGri_Dni_gas1)
        list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, self.n, amV, bmV, z_gas)
        list_PHI = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)
        ##############
        ##############
        # new
        Dew_P = self.calc_Dew_P(yi_list, list_PHI, list_gama)
        list_x = self.calc_x(yi_list, list_PHI, Dew_P, list_gama)

        amL = self.calc_am(EoS, list_x)
        bmL = self.calc_bm(EoS, list_x)
        z_liq = self.calc_z_liquid(EoS, amL, bmL, self.T, Dew_P)
        list_zi_liq = self.calc_z_i_liq(EoS, self.T, Dew_P)
        Gr_liq = self.calc_Gr(EoS, z_liq)
        list_Gr_i_liq = self.calc_Gr_i(EoS, list_zi_liq, self.T, Dew_P)
        list_DGri_Dni_liq = self.calc_DGri_Dni(amL, bmL, z_liq, list_x, EoS)
        list_GbarE_i = self.calc_GbarE_i(Gr_liq, list_Gr_i_liq, list_DGri_Dni_liq)
        list_gama = self.calc_gama(EoS, list_GbarE_i)

        Dew_P_new = self.calc_Dew_P(yi_list, list_PHI, list_gama)
        list_x_new = self.calc_x(yi_list, list_PHI, Dew_P, list_gama)

        amV = self.calc_am(EoS, yi_list)
        bmV = self.calc_bm(EoS, yi_list)
        z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, self.P)
        Gr_gas = self.calc_Gr(EoS, z_gas)
        # list_zi_gas = self.calc_z_i_gas(EoS, self.T, bubble_P)
        # list_Gr_i_gas = self.calc_Gr_i(EoS, list_zi_gas, self.T, bubble_P)
        list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, yi_list, EoS)
        list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, self.n, list_DGri_Dni_gas1)
        list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, self.n, amV, bmV, z_gas)
        list_PHI = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)

        #########
        #########
        while abs(Dew_P_new - Dew_P_old) > bubble_dew_pressure_acc:
            Dew_P_old = Dew_P_new
            Dew_P_new = self.calc_Dew_P(yi_list, list_PHI, list_gama)
            list_x = self.calc_x(yi_list, list_PHI, Dew_P_new, list_gama)

            amL = self.calc_am(EoS, list_x)
            bmL = self.calc_bm(EoS, list_x)
            z_liq = self.calc_z_liquid(EoS, amL, bmL, self.T, Dew_P_new)
            list_zi_liq = self.calc_z_i_liq(EoS, self.T, Dew_P_new)
            Gr_liq = self.calc_Gr(EoS, z_liq)
            list_Gr_i_liq = self.calc_Gr_i(EoS, list_zi_liq, self.T, Dew_P_new)
            list_DGri_Dni_liq = self.calc_DGri_Dni(amL, bmL, z_liq, list_x, EoS)
            list_GbarE_i = self.calc_GbarE_i(Gr_liq, list_Gr_i_liq, list_DGri_Dni_liq)
            list_gama = self.calc_gama(EoS, list_GbarE_i)

            Dew_P_new = self.calc_Dew_P(yi_list, list_PHI, list_gama)
            list_x = self.calc_x(yi_list, list_PHI, Dew_P_new, list_gama)

            amV = self.calc_am(EoS, yi_list)
            bmV = self.calc_bm(EoS, yi_list)
            z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, self.P)
            Gr_gas = self.calc_Gr(EoS, z_gas)
            # list_zi_gas = self.calc_z_i_gas(EoS, self.T, bubble_P)
            # list_Gr_i_gas = self.calc_Gr_i(EoS, list_zi_gas, self.T, bubble_P)
            list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, yi_list, EoS)
            list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, self.n, list_DGri_Dni_gas1)
            list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, self.n, amV, bmV, z_gas)
            list_PHI = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)

        return [Dew_P_new, list_PHI, list_gama, list_x]

    def calc_initial_PHI_list(self, dew_P, bubble_P, bubble_PHI_list, dew_PHI_list):
        cons = (self.P - dew_P) / (bubble_P - dew_P)
        new_PHI_list = []
        for i in range(len(bubble_PHI_list)):
            PHI_bub = bubble_PHI_list[i]
            PHI_dew = dew_PHI_list[i]
            new_PHI = (PHI_bub - PHI_dew) * cons + PHI_dew
            new_PHI_list.append(new_PHI)
        return new_PHI_list

    def calc_initial_gama_list(self, dew_P, bubble_P, bubble_gama_list, dew_gama_list):
        cons = (self.P - dew_P) / (bubble_P - dew_P)
        new_gama_list = []
        for i in range(len(bubble_gama_list)):
            gama_bub = bubble_gama_list[i]
            gama_dew = dew_gama_list[i]
            new_gama = (gama_bub - gama_dew) * cons + gama_dew
            new_gama_list.append(new_gama)
        return new_gama_list

    def calc_ki(self, gama_list, PHI_list):
        ki_list = []
        for i in range(len(gama_list)):
            component = self.components[i]
            Psat = component.Psat
            gama = gama_list[i]
            PHI = PHI_list[i]
            ki = gama * Psat / (PHI * self.P)
            ki_list.append(ki)
        return ki_list

    def calc_F_flash(self, ki_list, V):
        F = 0
        for i in range(len(ki_list)):
            zi = self.compositions[i]
            ki = ki_list[i]
            Fi = zi * (ki - 1) / (1 + V * (ki - 1))
            F += Fi
        return F

    def calc_DF_DV_flash(self, ki_list, V):
        DF_DV = 0
        for i in range(len(ki_list)):
            zi = self.compositions[i]
            ki = ki_list[i]
            Fi = zi * ((ki - 1) ** 2) / ((1 + V * (ki - 1)) ** 2)
            DF_DV += Fi
        return - DF_DV

    def newtone_V_calc(self, ki_list, initial_V):
        F_old = self.calc_F_flash(ki_list, initial_V)
        DF_DV_old = self.calc_DF_DV_flash(ki_list, initial_V)
        V_old = initial_V

        V_new = V_old - (F_old / DF_DV_old)
        F_new = self.calc_F_flash(ki_list, V_new)
        DF_DV_new = self.calc_DF_DV_flash(ki_list, V_new)

        acc = 10 ** - 4
        while abs(V_old - V_new) > acc:
            V_old = V_new
            F_old = F_new
            DF_DV_old = DF_DV_new
            V_new = V_old - (F_old / DF_DV_old)
            F_new = self.calc_F_flash(ki_list, V_new)
            DF_DV_new = self.calc_DF_DV_flash(ki_list, V_new)
        return V_new

    def calc_flash_xi_yi(self, V, ki_list):
        xi_list = []
        yi_list = []
        for i in range(len(ki_list)):
            zi = self.compositions[i]
            ki = ki_list[i]
            xi = zi / (1 + V * (ki - 1))
            yi = ki * xi
            xi_list.append(xi)
            yi_list.append(yi)
        return xi_list, yi_list

    def flash_calc(self, EoS, bubble_dew_pressure_acc):
        phase_list = self.phase_list
        bubble_props = self.calc_bubble_point(EoS, bubble_dew_pressure_acc)
        bubble_P = bubble_props[0]
        bubble_PHI_list = bubble_props[1]
        bubble_gama_list = bubble_props[2]

        dew_props = self.calc_dew_point(EoS, bubble_dew_pressure_acc)
        dew_P = dew_props[0]
        dew_PHI_list = dew_props[1]
        dew_gama_list = dew_props[2]

        if self.P > bubble_P:
            liquid_phase = LiquidPhase(self, self.compositions, EoS, self.n)
            self.__phase_list.clear()
            phase_list.append(liquid_phase)
            # liquid phase must be created
            # composition of liquid phase is self.composition
        elif self.P < dew_P:
            gas_phase = GasPhase(self, self.compositions, EoS, self.n)
            self.__phase_list.clear()
            phase_list.append(gas_phase)
            # gas phase must be created
            # composition = self.composition
        else:
            # we have liquid and gas phase
            # here we calc V, xi, yi
            initial_PHI_list = self.calc_initial_PHI_list(dew_P, bubble_P, bubble_PHI_list, dew_PHI_list)
            initial_gama_list = self.calc_initial_gama_list(dew_P, bubble_P, bubble_gama_list, dew_gama_list)
            initial_V = (bubble_P - self.P) / (bubble_P - dew_P)
            ki_list = self.calc_ki(initial_gama_list, initial_PHI_list)

            V = self.newtone_V_calc(ki_list, initial_V)
            old_V = initial_V
            xi_list, yi_list = self.calc_flash_xi_yi(V, ki_list)
            acc = 10 ** -4
            while abs(V - old_V) > acc:

                amL = self.calc_am(EoS, xi_list)
                bmL = self.calc_bm(EoS, xi_list)
                z_liq = self.calc_z_liquid(EoS, amL, bmL, self.T, self.P)
                list_zi_liq = self.calc_z_i_liq(EoS, self.T, self.P)
                Gr_liq = self.calc_Gr(EoS, z_liq)
                list_Gr_i_liq = self.calc_Gr_i(EoS, list_zi_liq, self.T, self.P)
                list_DGri_Dni_liq = self.calc_DGri_Dni(amL, bmL, z_liq, xi_list, EoS)
                list_GbarE_i = self.calc_GbarE_i(Gr_liq, list_Gr_i_liq, list_DGri_Dni_liq)

                list_gama = self.calc_gama(EoS, list_GbarE_i)

                amV = self.calc_am(EoS, yi_list)
                bmV = self.calc_bm(EoS, yi_list)
                z_gas = self.calc_z_gas(EoS, amV, bmV, self.T, self.P)
                Gr_gas = self.calc_Gr(EoS, z_gas)
                # list_zi_gas = self.calc_z_i_gas(EoS, self.T, bubble_P)
                # list_Gr_i_gas = self.calc_Gr_i(EoS, list_zi_gas, self.T, bubble_P)
                list_DGri_Dni_gas1 = self.calc_DGri_Dni(amV, bmV, z_gas, yi_list, EoS)
                list_phi_i_hat = self.calc_phi_i_hat(Gr_gas, self.T, EoS, (V * self.n), list_DGri_Dni_gas1)
                list_phi_i_hat_sat = self.calc_ln_phi_i_hat_sat(Gr_gas, self.T, EoS, (V * self.n), amV, bmV, z_gas)

                list_PHI = self.calc_PHI(list_phi_i_hat, list_phi_i_hat_sat)

                ki_list = self.calc_ki(list_gama, list_PHI)
                old_V = V
                V = self.newtone_V_calc(ki_list, initial_V)
                xi_list, yi_list = self.calc_flash_xi_yi(V, ki_list)
            gas_phase = GasPhase(self, yi_list, EoS, V * self.n)
            liquid_phase = LiquidPhase(self, xi_list, EoS, (1 - V) * self.n)
            phase_list.append(gas_phase)
            phase_list.append(liquid_phase)
        return phase_list
