import pickle

import Component


class PureComponent(Component):
    all_props = dict()
    all_pure_components = []

    @staticmethod
    def load_component(component_name):
        with open("components\\" + component_name + ".dat", "rb") as f:
            return pickle.load(f)

    def __init__(self, name, Tc, Pc):
        self.__name = name
        self.__tc = Tc
        self.__Pc = Pc
        PureComponent.all_pure_components.append(self)

    # in this class we use serialization for some very usable components such as water etc.
    # at the end of this method the component must be added to all pure components list
    def save(self):
        with open("components\\" + self.name, "wb") as file:
            pickle.dump(self, file)


methane = PureComponent("methane", 200, 300)
ethane = PureComponent("ethane", 500, 700)

methane.save()
ethane.save()
