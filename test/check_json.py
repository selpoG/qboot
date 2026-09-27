"""Check PMP JSON independently of the C++ writer, using only Python's stdlib."""

import json
import os
from decimal import Decimal, getcontext
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


def close(actual, expected):
    assert abs(actual - expected) < Decimal("1e-60"), (actual, expected)


def read_matrix(tokens):
    rows, columns = int(next(tokens)), int(next(tokens))
    return [[Decimal(next(tokens)) for _ in range(columns)] for _ in range(rows)]


def consistency(directory, degree):
    name = f"consistency-{degree}"
    pmp = json.loads((directory / f"{name}.json").read_text())
    xml = ET.parse(directory / f"{name}.xml").getroot()
    direct = directory / f"{name}-sdp"
    objective = numbers(pmp["objective"])
    for actual, expected in zip(objective, [Decimal("11.625"), Decimal(5)]):
        close(actual, expected)
    assert objective == [Decimal(x.text) for x in xml.find("objective")]
    raw = (direct / "objectives").read_text().split()
    assert int(raw[1]) == 1
    close(Decimal(raw[0]), objective[0])
    close(Decimal(raw[2]), objective[1])
    block = pmp["PositiveMatrixWithPrefactorArray"][0]
    xml_block = xml.find("polynomialVectorMatrices")[0]
    polynomials = numbers(block["polynomials"])
    xml_polynomials = iter(xml_block.find("elements"))
    for r, row in enumerate(polynomials):
        for c, vector in enumerate(row):
            factor = r + 2 if r == c else 1
            reference = next(xml_polynomials)
            for n, polynomial in enumerate(vector):
                assert len(polynomial) == degree + 1
                for j, value in enumerate(polynomial):
                    expected = Decimal("0.5") + Decimal("1.5") * j if n == 0 else Decimal(2 + j)
                    close(value, factor * expected)
                    close(value, Decimal(reference[n][j].text))
    points = numbers(block["samplePoints"])
    scales = numbers(block["sampleScalings"])
    assert points == list(range(1, degree + 2))
    assert scales == list(range(2, degree + 3))
    for field in ("samplePoints", "sampleScalings"):
        assert numbers(block[field]) == [Decimal(x.text) for x in xml_block.find(field)]
    b = read_matrix(iter((direct / "free_var_matrix.0").read_text().split()))
    raw_c = (direct / "primal_objective_c.0").read_text().split()
    assert int(raw_c[0]) == 3 * (degree + 1)
    rhs = iter(map(Decimal, raw_c[1:]))
    rows = iter(b)
    for r in range(2):
        for c in range(r + 1):
            for x, scale in zip(points, scales):
                const, variable = [sum(coefficient * x**j for j, coefficient in enumerate(p))
                                   for p in polynomials[r][c]]
                close(next(rhs), const * scale)
                close(next(rows)[0], -variable * scale)
    tokens = iter((direct / "bilinear_bases.0").read_text().split())
    assert int(next(tokens)) == 1
    for parity in (0, 1):
        basis = read_matrix(tokens)
        expected_rows = degree // 2 + 1 if parity == 0 else (degree + 1) // 2
        assert len(basis) == expected_rows, (degree, parity, len(basis), expected_rows)
        assert len(block[f"bilinearBasis_{parity}"]) == expected_rows
        for m, row in enumerate(basis):
            for value, x, scale in zip(row, points, scales):
                close(value * value, x ** (2 * m + parity) * scale)
    assert next(tokens, None) is None


def check(directory):
    getcontext().prec = 100
    for degree in range(3):
        consistency(directory, degree)
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
    # and the two constraints are 3(z-3)/2 >= 0 and 5-z >= 0.
    assert numbers(pmp["objective"]) == [Decimal("-0.375"), Decimal("0.5")]
    assert numbers(json_blocks[0]["polynomials"])[0][0] == [[Decimal("-4.5")], [Decimal("1.5")]]
    assert numbers(json_blocks[1]["polynomials"])[0][0] == [[5], [-1]]


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="qboot-json-") as temporary:
        directory = Path(temporary)
        subprocess.run([sys.argv[1]], env=dict(os.environ, QBOOT_TEST_OUTPUT=str(directory)), check=True)
        check(directory)
