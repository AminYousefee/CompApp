class PropertyContainer:
    required_property = dict()

    def __init__(self, properties):
        self.properties = properties

    def get_values(self, properties):
        pass

    # properties is a list of name of properties

    def set_values(self, properrties):
        pass

    # properties is a dict of properties which key is prop name

    def calculate_prop(self, prop_name):
        pass

    def update_prop(self, prop_name):
        pass

    def get_prop_names(self):
        pass
    # returns keys of required property (names)
