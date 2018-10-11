class rocket:

    def __init__(self, body, fin, nosecone, motor, recovery, payload, slug):
        self.body = body
        self.fin = fin
        self.nosecone = nosecone
        self.motor = motor
        self.recovery = recovery
        self.payload = payload
        self.slug = slug

class body:

    def __init__(self, d, l):
        self.d = d
        self.l = l

    @staticmethod
    def from_file(file):
        d = find_parameter(file, "diameter")
        l = find_parameter(file, "length")
        return body(d, l)


def find_parameter(file, parameter):
    file = open(file, 'r')
    arr = ["", ""]
    while arr[0] != parameter.lower():
        base = file.readline()
        base = base.replace(" ", "")
        arr = base.split("=")
    file.close()
    return eval(arr[1])

