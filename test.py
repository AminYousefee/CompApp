class Test:
    counter = 0

    def __init__(self, a, b):
        self.a = a
        self.b = b
        Test.counter += 1

    @staticmethod
    def get_count():
        return Test.counter


a = Test(1, 2)
b = Test(2, 3)

print(Test.counter)
