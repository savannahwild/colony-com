
#Y = Species("Y",U_Y)
class Species:
    def __init__(self, name, U):
        self.name = name
        #Y.name = "Y"
        self.U = U
        #Y.U = U_Y
        self.behaviour = None
    
    def get_U(self):
        return self.U
        #Y.get_U returns U_Y

    #Y.set_U(U)
    def set_U(self, U):
        self.U = U
        #returns Y.U = U

    def get_name(self):
        return self.name
        #Y.get_name returns Y.name which is "Y"
    
    #Y.set_behaviour(Y_behaviour)
    def set_behaviour(self, behaviour):
        self.behaviour = behaviour
        #Y.set_behaviour = Y_behaviour
        