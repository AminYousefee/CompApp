from Just_name import Just_name


class Name_Last(Just_name):
    counter = 0

    def say_hello(self):
        print("hello" + self.name + " " + self.last_name)

    def __init__(self, name, last_name):
        super().__init__(name)
        self.last_name = last_name


"""
    @property
    def last_name(self):
        print("getter of name")
        return self.__last_name

    @last_name.setter
    def last_name(self, val):
        self.__last_name = val
        print("name setter")
"""
