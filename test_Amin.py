from EoS import EoS
from Fluid import Fluid
from PureComponent import PureComponent

PR = EoS.create_pr()

hexane = {
    "name": "hexane",
    "Tc": 507.6,
    "Pc": 30.25 * (10 ** 5),
    "omega": 0.301,
    "Tnb": 341.9,
    "MW": 86.177,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * 3.025, EoS.R * 53.722 * 10 ** -3, EoS.R * (-16.791 * 10 ** -6), EoS.R * 0],
    "Antoine coeff": [13.8193, 2696.04, 224.317]
}

water = {
    "name": "water",
    "Tc": 647.1,
    "Pc": 220.55 * (10 ** 5),
    "omega": 0.345,
    "Tnb": 373.2,
    "MW": 18.015,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * 3.47, EoS.R * 1.45 * 10 ** -3, EoS.R * 0, EoS.R * (0.121 * (10 ** 5))],
    "Antoine coeff": [16.3872, 3885.7, 230.17]
}

ethanol = {
    "name": "ethanol",
    "Tc": 513.9,
    "Pc": 61.48 * (10 ** 5),
    "omega": 0.645,
    "Tnb": 351.4,
    "MW": 46.069,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * 3.518, EoS.R * 20.001 * 10 ** -3, EoS.R * (-6.002 * 10 ** -6), EoS.R * 0],
    "Antoine coeff": [16.8958, 3795.17, 230.918]
}

methanol = {
    "name": "methanol",
    "Tc": 512.6,
    "Pc": 80.97 * (10 ** 5),
    "omega": 0.564,
    "Tnb": 337.9,
    "MW": 32.042,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * 2.211, EoS.R * 12.216 * 10 ** -3, EoS.R * (-3.45 * 10 ** -6), EoS.R * 0],
    "Antoine coeff": [16.5785, 3638.27, 239.5]
}

acetone = {
    "name": "acetone",
    "Tc": 508.2,
    "Pc": 47.01 * (10 ** 5),
    "omega": 0.307,
    "Tnb": 329.4,
    "MW": 58.08,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * 3.025, EoS.R * 53.722 * 10 ** -3, EoS.R * (-16.791 * 10 ** -6), EoS.R * 0],
    "Antoine coeff": [14.3145, 2756.22, 228.06]
}

benzene = {
    "name": "benzene",
    "Tc": 562.2,
    "Pc": 48.98 * (10 ** 5),
    "omega": 0.210,
    "Tnb": 353.2,
    "MW": 78.114,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * (-0.206), EoS.R * 39.064 * 10 ** -3, EoS.R * (-13.301 * 10 ** -6), EoS.R * 0],
    "Antoine coeff": [13.7819, 2726.81, 217.572]
}

toluene = {
    "name": "toluene",
    "Tc": 591.8,
    "Pc": 41.06 * (10 ** 5),
    "omega": 0.262,
    "Tnb": 383.8,
    "MW": 92.141,
    "viscosity25": 0,
    "Cp coefficients": [EoS.R * 2.050, EoS.R * 12.216 * 10 ** -3, EoS.R * (-3.45 * 10 ** -6), EoS.R * 0],
    "Antoine coeff": [13.9320, 3056.96, 217.625]
}

Toluene = PureComponent(toluene)
# PureComponent.save_component(toluene)
# Toluene1 = PureComponent.load_component("toluene")

Benzene = PureComponent(benzene)
# PureComponent.save_component(benzene)
# Benzene1 = PureComponent.load_component("benzene")

Water = PureComponent(water)
# PureComponent.save_component(water)
# Water1 = PureComponent.load_component("water")

fluid1 = Fluid(2, 288, 101325, [Benzene, Toluene], [0.6, 0.4], PR)
fluid2 = Fluid(2, 363, 101325, [Water], [1], PR)

phase_list1 = fluid1.phase_list
phase_list2 = fluid2.phase_list

# here you can read all phase properties like: phase_list1[0].H
new_conditions1 = {"T": 300, "P": 100000, "composition": [0.5, 0.5]}
fluid1.update(new_conditions1)

new_conditions2 = {"T": 300, "P": 100000}
fluid2.update(new_conditions2)

new_conditions3 = {"T": 320}
fluid1.update(new_conditions3)

phase_list1_2 = fluid1.phase_list
phase_list2_2 = fluid2.phase_list
