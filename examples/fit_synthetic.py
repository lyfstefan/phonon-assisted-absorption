"""Synthetic fit demonstration, not a reproduction of experimental measurements."""

from pathlib import Path
from phonon_absorption.demo import run

if __name__ == "__main__":
    print(run(Path("examples/output")))
