import Component
import pickle


class PureComponent(Component):
    all_props = dict()
    all_pure_components = []

    @staticmethod
    def load_component(component_name):
        with open("components\\" + component_name + ".dat", "rb") as f:
            return pickle.load(f)

    def __init__(self):
        pass

    # in this class we use serialization for some very usable components such as water etc.
    # at the end of this method the component must be added to all pure components list
    def __save(self):
        with open("components\\" + self.name, "wb") as file:
            pickle.dump(self, file)
