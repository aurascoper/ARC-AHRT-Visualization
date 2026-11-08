"""ESTAR table: verified NIST electron stopping-power / CSDA-range data, hash-pinned.
Same provenance discipline as XCOM. Log-log interpolation (standard for stopping tables).
"""
from dataclasses import dataclass
import csv, math, os

@dataclass
class EstarRow:
    energy_MeV: float
    collision: float   # MeV cm2/g
    radiative: float
    total: float
    csda_g_cm2: float | None = None  # only on the range grid

class EstarTable:
    def __init__(self, stop_rows, range_rows):
        self.stop = sorted(stop_rows, key=lambda r: r.energy_MeV)
        self.rng = sorted(range_rows, key=lambda r: r.energy_MeV)

    @staticmethod
    def _log_interp(E, xs, ys):
        lx = math.log(E)
        for i in range(1, len(xs)):
            if xs[i] >= E:
                t = (lx - math.log(xs[i-1])) / (math.log(xs[i]) - math.log(xs[i-1]))
                return math.exp(math.log(ys[i-1]) + t * (math.log(ys[i]) - math.log(ys[i-1])))
        raise ValueError(f"energy {E} MeV outside ESTAR grid [{xs[0]}, {xs[-1]}]")

    def stopping(self, E_MeV):
        xs = [r.energy_MeV for r in self.stop]
        return (self._log_interp(E_MeV, xs, [r.collision for r in self.stop]),
                self._log_interp(E_MeV, xs, [r.radiative for r in self.stop]),
                self._log_interp(E_MeV, xs, [r.total for r in self.stop]))

    def csda_range_cm(self, E_MeV, density_g_cm3=1.0):
        xs = [r.energy_MeV for r in self.rng]
        g_cm2 = self._log_interp(E_MeV, xs, [r.csda_g_cm2 for r in self.rng])
        return g_cm2 / density_g_cm3, g_cm2

def load_estar(stop_csv, range_csv):
    def _rows(path, csda):
        rows = []
        with open(path) as f:
            for ln in f:
                if ln.startswith("#") or not ln.strip():
                    continue
                v = [float(x) for x in ln.strip().split(",")]
                rows.append(EstarRow(v[0], v[1], v[2], v[3], v[4] if csda else None))
        return rows
    return EstarTable(_rows(stop_csv, False), _rows(range_csv, True))
