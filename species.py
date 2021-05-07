#Species class
#Define behaviour

class Species:
    def __init__(self, name, U):
        self.name = name
        self.U = U
        self.behaviour = None

    def get_U(self):
        return self.U

    def set_U(self, U):
        self.U = U

    def get_name(self):
        return self.name

    def set_behaviour(self, behaviour):
        self.behaviour = behaviour
