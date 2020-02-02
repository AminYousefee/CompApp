from Phase import Phase


class GasPhase(Phase):
    def __init__(self, fluid, compositions, EoS):
        super().__init__(fluid, compositions, EoS)
