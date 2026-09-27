"""Check a solved qboot PMP against its analytic optimum."""

from decimal import Decimal
from pathlib import Path
import re
import sys

output = Path(sys.argv[1]).read_text()
assert '"found primal-dual optimal solution"' in output, output
for field in ("primalObjective", "dualObjective"):
    match = re.search(rf"{field}\s*=\s*([^;]+);", output)
    assert match, output
    assert abs(Decimal(match[1]) - Decimal("2.125")) < Decimal("1e-25"), output
