from Phase import Phase


class LiquidPhase(Phase):
    def __init__(self, fluid, compositions, EoS):
        super().__init__(fluid, compositions, EoS)
