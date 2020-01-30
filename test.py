class Test:
    counter = 0

    def __init__(self, name):
        self.__name = name

    @property
    def name(self):
        print("getter of name")
        return self.__name

    @name.setter
    def name(self, val):
        print("name setter")


a = Test("Mina")
print(a.name)
a.name = "Amin"
print(a.name)
