import Component


class PureComponent(Component):
    all_props = dict()

    @staticmethod
    def load_component(component_name):
        with open(component_name + ".dat", "rb"):

    def __init__(self):
        pass
    # in this class we use serialization for some very usable components such as water etc.
