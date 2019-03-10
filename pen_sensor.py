class VirtualSensor:
    def __init__(self):
        self.heigth = 0
        self.acceleration = 0
        pass

    def in_heigth(self, heigth):
        self.height = heigth

    def in_acceleration(self, acceleration):
        self.acceleration = acceleration

    def get_height(self):
        return self.heigth

    def get_acceleration(self):
        return self.acceleration
