import pickle


class PureComponent:
    all_props = dict()
    all_pure_components = dict()

    @staticmethod
    def load_components():
        with open("components\\components.dat", "rb") as f:
            a = pickle.load(f)
            PureComponent.all_pure_components = a
            return a

    def __init__(self, name, Tc, Pc):
        self.name = name
        self.tc = Tc
        self.Pc = Pc
        PureComponent.all_pure_components[self.name] = self

    # in this class we use serialization for some very usable components such as water etc.
    # at the end of this method the component must be added to all pure components list
    @staticmethod
    def save():
        with open("components\\components.dat", "wb") as file:
            pickle.dump(PureComponent.all_pure_components, file)


methane = PureComponent("methane", 200, 300)
ethane = PureComponent("ethane", 500, 700)
PureComponent.save()
PureComponent.all_pure_components = None
all = PureComponent.load_components()
print(all)
print(all["ethane"])

if PureComponent.all_pure_components != None:
    print("ye machi az erfan bigir")
