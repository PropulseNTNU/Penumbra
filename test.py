class parent:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def func1(self):
        print(self.x, self.y)

class child(parent):
    def __init__(self, x, y):
        parent.__init__(self, x, y)

    def func2(self):
        print(self.x, self.y)

childObj = child(10, 20)
childObj.func2()
