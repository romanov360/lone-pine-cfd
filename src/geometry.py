"""Bottle and environment definitions."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from props import GLASS, Solid, water


@dataclass(frozen=True)
class Bottle:
    """A cylindrical bottle: liquid core inside a solid shell."""

    D_outer: float = 0.070       # m, outside diameter
    H_outer: float = 0.145       # m, outside height
    wall: float = 0.0030         # m, wall thickness
    shell: Solid = GLASS
    name: str = "500 mL glass bottle"

    @property
    def D_inner(self) -> float:
        return self.D_outer - 2.0 * self.wall

    @property
    def H_inner(self) -> float:
        return self.H_outer - 2.0 * self.wall

    @property
    def V_liquid(self) -> float:
        return math.pi * (self.D_inner / 2.0) ** 2 * self.H_inner

    @property
    def V_shell(self) -> float:
        outer = math.pi * (self.D_outer / 2.0) ** 2 * self.H_outer
        return outer - self.V_liquid

    @property
    def A_outer(self) -> float:
        """Total wetted outside area (side + both ends)."""
        side = math.pi * self.D_outer * self.H_outer
        ends = 2.0 * math.pi * (self.D_outer / 2.0) ** 2
        return side + ends

    @property
    def A_inner(self) -> float:
        side = math.pi * self.D_inner * self.H_inner
        ends = 2.0 * math.pi * (self.D_inner / 2.0) ** 2
        return side + ends

    @property
    def A_side_outer(self) -> float:
        return math.pi * self.D_outer * self.H_outer

    def C_liquid(self, T_C: float = 25.0) -> float:
        """Heat capacity of the contained water, J/K."""
        w = water(T_C)
        return self.V_liquid * w.rho * w.cp

    def C_shell(self) -> float:
        return self.V_shell * self.shell.rho * self.shell.cp

    def C_total(self, T_C: float = 25.0) -> float:
        return self.C_liquid(T_C) + self.C_shell()

    @property
    def R_wall(self) -> float:
        """
        Conduction resistance of the shell, K/W, using the exact cylindrical
        form for the side wall in parallel with plane slabs for the ends.
        """
        r_i, r_o = self.D_inner / 2.0, self.D_outer / 2.0
        k = self.shell.k
        R_side = math.log(r_o / r_i) / (2.0 * math.pi * k * self.H_inner)
        A_end = math.pi * r_i**2
        R_ends = self.wall / (k * A_end) / 2.0     # two ends in parallel
        return 1.0 / (1.0 / R_side + 1.0 / R_ends)

    @property
    def U_wall(self) -> float:
        """Shell conductance expressed as an equivalent film coefficient on
        the outside area, W/m^2-K. Lets it be compared with h directly."""
        return 1.0 / (self.R_wall * self.A_outer)

    def describe(self) -> str:
        return (
            f"{self.name}\n"
            f"  outside      {self.D_outer*1000:.0f} mm x {self.H_outer*1000:.0f} mm\n"
            f"  wall         {self.wall*1000:.1f} mm {self.shell.name}"
            f"  (k={self.shell.k} W/m-K)\n"
            f"  liquid       {self.V_liquid*1e6:.0f} mL\n"
            f"  wetted area  {self.A_outer*1e4:.0f} cm^2\n"
            f"  C_liquid     {self.C_liquid()/1000:.2f} kJ/K"
            f"   C_shell {self.C_shell()/1000:.2f} kJ/K"
            f"   ({self.C_shell()/self.C_total()*100:.0f}% in the glass)\n"
            f"  R_wall       {self.R_wall*1000:.3f} mK/W"
            f"  -> equivalent h {self.U_wall:.0f} W/m^2-K"
        )


# The reference bottle sized so the liquid volume lands on 500 mL.
def _sized_500ml() -> Bottle:
    b = Bottle()
    target = 500e-6
    lo, hi = 0.10, 0.25
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if replace(b, H_outer=mid).V_liquid < target:
            lo = mid
        else:
            hi = mid
    return replace(b, H_outer=0.5 * (lo + hi))


BOTTLE_500 = _sized_500ml()


@dataclass(frozen=True)
class StillBox:
    """A tub of motionless cold water."""

    volume: float = 0.020       # m^3 (20 L)
    T0: float = 10.0            # C
    insulated: bool = True      # no loss to the room
    name: str = "box"

    @property
    def side(self) -> float:
        """Cube-equivalent side length, m."""
        return self.volume ** (1.0 / 3.0)

    def C(self) -> float:
        w = water(self.T0)
        return self.volume * w.rho * w.cp


@dataclass(frozen=True)
class Creek:
    """Flowing water: effectively infinite, at a fixed temperature."""

    U: float = 0.35             # m/s, depth-averaged approach velocity
    T0: float = 10.0            # C
    void_fraction: float = 0.0  # entrained air by volume
    turbulence_intensity: float = 0.08
    name: str = "creek"

    def C(self) -> float:
        return math.inf


if __name__ == "__main__":
    print(BOTTLE_500.describe())
    print()
    for V in (0.005, 0.010, 0.020, 0.060, 0.200, 1.0):
        box = StillBox(volume=V)
        print(f"box {V*1000:7.1f} L : C={box.C()/1000:8.1f} kJ/K"
              f"   C_box/C_bottle = {box.C()/BOTTLE_500.C_total():6.1f}"
              f"   side={box.side*100:5.1f} cm")
