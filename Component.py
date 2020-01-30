import ProperyContainer


class Component(ProperyContainer):
    required_props_name = ["Tc", "Pc", "omega", "MW", "Tnb"]
    all_props = dict()
    """
    Tc
    Pc
    omega
    MW
    viscosity at 25
    nbp
    """
    """
    a and b
    a , b , c , d must be given to us
    thermal conductivity
    """

    def __init__(self, initial_prop_dict):
        if len(Component.required_props_name) != len(initial_prop_dict):
            raise
        for prop_name in Component.required_props_name:
            if prop_name in initial_prop_dict:
                if prop_name == "Tc":
                    self.Tc = initial_prop_dict["Tc"]
                if prop_name == "Pc":
                    self.Pc = initial_prop_dict["Pc"]
                if prop_name == "omega":
                    self.omega = initial_prop_dict["omega"]
                if prop_name == "Tnb":
                    self.Tnb = initial_prop_dict["Tnb"]
                if prop_name == "MW":
                    self.MW = initial_prop_dict["MW"]
                if prop_name == "viscosity25":
                    self.visco25 = initial_prop_dict["viscosity25"]
                if prop_name == "Cp coefficients":
                    self.Cp_coeff = initial_prop_dict["Cp coefficients"]

        # now other properties must be calculate from initial properties
