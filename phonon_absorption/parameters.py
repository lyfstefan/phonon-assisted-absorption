"""Rounded Table I / thesis Table 4 anchors; not new fits or validated domains."""
from .models import Channel

ANCHORS={
    "I":dict(phonon_energy_eV=.0346,linewidth_eV=.0241,kBT_eV=.0284,
             channels=[Channel(1.281,1.5755,.3728**2),Channel(.897,2.6136,.4692**2)]),
    "II":dict(phonon_energy_eV=.0320,linewidth_eV=.0231,kBT_eV=.0232,
              channels=[Channel(1.281,.85,.45**2,.60),Channel(.898,1.34,.51**2,.48)])
}
