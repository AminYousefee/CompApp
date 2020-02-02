from Phase import Phase


class GasPhase(Phase):
    def __init__(self, fluid, compositions):
        super().__init__(fluid)
        self.__compositions = compositions
