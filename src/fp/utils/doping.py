#region modules
#endregions

#region variables
#endregions

#region functions
#endregions

#region classes
class Doping:
    def __init__(
        self,
        mixing_atoms: list,
        replacement_atom: str,
        fraction: float,
        **kwargs,
    ):
        self.mixing_atoms: list = mixing_atoms
        self.replacement_atom: str = replacement_atom
        self.fraction: float = fraction

        for key, value in kwargs.items():
            setattr(self, key, value)

    def gen_pseudo(self):
        pass
#endregions
