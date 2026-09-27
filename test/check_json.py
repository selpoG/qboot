"""Check PMP JSON independently of the C++ writer, using only Python's stdlib."""

import json
from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def numbers(value):
    if isinstance(value, list):
        return [numbers(item) for item in value]
    assert isinstance(value, str), "SDPB numbers must be strings"
    result = Decimal(value)
    assert result.is_finite()
    return result


def check(directory):
    matrices = json.loads((directory / "matrices.json").read_text())
    assert numbers(matrices["objective"]) == [Decimal("0.125"), 1]
    assert numbers(matrices["normalization"]) == [1, 0]
    blocks = matrices["PositiveMatrixWithPrefactorArray"]
    assert len(blocks) == 3
    for degree, block in enumerate(blocks):
        assert numbers(block["samplePoints"]) == list(range(1, degree + 2))
        assert numbers(block["sampleScalings"]) == [1] * (degree + 1)
        assert len(block["bilinearBasis_0"]) == degree // 2 + 1
        assert len(block["bilinearBasis_1"]) == (degree + 1) // 2
        for parity in (0, 1):
            for k, polynomial in enumerate(numbers(block[f"bilinearBasis_{parity}"])):
                assert polynomial == [0] * k + [1]
        polynomials = numbers(block["polynomials"])
        assert len(polynomials) == 2
        for r, row in enumerate(polynomials):
            assert len(row) == 2
            for c, vector in enumerate(row):
                assert len(vector) == 2
                assert all(len(p) == degree + 1 for p in vector)
                assert all(x == 0 for p in vector for x in p[1:])
                if r == c:
                    assert abs(vector[0][0] - Decimal("1.23456789012345678901234567890123456789")) < Decimal("1e-70")
                    assert vector[1][0] == -1
                else:
                    assert all(x == 0 for p in vector for x in p)

    source = (directory / "bounded.json").read_text()
    assert source == (directory / "bounded-parallel.json").read_text()
    pmp = json.loads(source)
    xml = ET.parse(directory / "bounded.xml").getroot()
    assert numbers(pmp["objective"]) == [Decimal(x.text) for x in xml.find("objective")]
    json_blocks = pmp["PositiveMatrixWithPrefactorArray"]
    xml_blocks = xml.find("polynomialVectorMatrices")
    assert len(json_blocks) == len(xml_blocks) == 2
    for block, reference in zip(json_blocks, xml_blocks):
        for field in ("samplePoints", "sampleScalings"):
            assert numbers(block[field]) == [Decimal(x.text) for x in reference.find(field)]
        actual = [p for row in numbers(block["polynomials"]) for p in row]
        expected = [[[Decimal(x.text) for x in p] for p in v] for v in reference.find("elements")]
        assert actual == expected

    # After eliminating y using z = 2y + 1: objective = -3/8 + z/2,
    # and the two constraints are (z-1)/2 >= 0 and 5-z >= 0.
    assert numbers(pmp["objective"]) == [Decimal("-0.375"), Decimal("0.5")]
    assert numbers(json_blocks[0]["polynomials"])[0][0] == [[Decimal("-0.5")], [Decimal("0.5")]]
    assert numbers(json_blocks[1]["polynomials"])[0][0] == [[5], [-1]]


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="qboot-json-") as temporary:
        directory = Path(temporary)
        subprocess.run([sys.argv[1], str(directory)], check=True)
        check(directory)
