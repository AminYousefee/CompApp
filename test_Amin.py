import pickle
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
        self.__name = val
        print("name setter")


a = Test("Mina")
print(a.name)
a.name = "Amin"
print(a.name)
b = Test("mammad")
print(b.name)
with open("pickle_test.dat", "wb") as file:
    pickle.dump(a, file)
    print("a is pickled")

with open("pickle_test.dat", "rb") as f:
    b = pickle.load(f)
    print("b ia like a and it's name must be Amin")
    print(b.name)

print(b.name)
