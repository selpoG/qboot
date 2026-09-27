"""Compare qboot's direct SDP with SDPB's conversions of XML and JSON PMP."""

import json
from decimal import Decimal, getcontext
from pathlib import Path
import sys

from check_json import close, numbers, read_matrix


def compare(actual, expected):
    if isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected)
        for a, b in zip(actual, expected):
            compare(a, b)
    else:
        close(actual, expected)


def check(directory):
    getcontext().prec = 100
    for degree in range(3):
        direct = directory / f"consistency-{degree}-sdp"
        objective = (direct / "objectives").read_text().split()
        b = read_matrix(iter((direct / "free_var_matrix.0").read_text().split()))
        c = list(map(Decimal, (direct / "primal_objective_c.0").read_text().split()[1:]))
        tokens = iter((direct / "bilinear_bases.0").read_text().split())
        assert int(next(tokens)) == 1
        even, odd = read_matrix(tokens), read_matrix(tokens)
        for format in ("json", "xml"):
            converted = directory / f"converted-{degree}-{format}"
            objectives = json.loads((converted / "objectives.json").read_text())
            close(numbers(objectives["constant"]), Decimal(objective[0]))
            compare(numbers(objectives["b"]), list(map(Decimal, objective[2:])))
            block = json.loads((converted / "block_data_0.json").read_text())
            for key, expected in (("B", b), ("c", c), ("bilinear_bases_even", even), ("bilinear_bases_odd", odd)):
                compare(numbers(block[key]), expected)


if __name__ == "__main__":
    check(Path(sys.argv[1]))
