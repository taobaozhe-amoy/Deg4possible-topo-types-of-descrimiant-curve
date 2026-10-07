# Deg4possible-topo-types-of-descrimiant-curve



To check 6 examples of topological types of deg 4 descrimiant curve of our spefici model.
Please run the following in the sage kenel while keep "quatic_fixed.sage" under the same folder
from pathlib import Path
import re
import builtins

# ============================================================
# ============================================================

VINZANT_FILE = Path("quartic_fixed.sage")

if not VINZANT_FILE.exists():
    raise FileNotFoundError(
        ""
        ""
    )

vinzant_source = VINZANT_FILE.read_text()


# ============================================================
# 
# ============================================================

quartics = [

    (
        "1. Empty set",
        """
(2*x^2 + 3*y^2 + 4*z^2)^2
- (x^2 + y^2 - z^2)
  *((-9*x^2 - 8*y^2 + 13*z^2)/10)
"""
    ),

    (
        "2. One oval",
        """
(z^2/10)^2
- (x^2 + y^2 - z^2)
  *(-x^2 - 2*y^2 - z^2)
"""
    ),

    (
        "3. Two non-nested ovals",
        """
(3*x*z)^2
- (x^2 + y^2 - z^2)
  *(-x^2 - 2*y^2 + 5*z^2)
"""
    ),

    (
        "4. Two nested ovals",
        """
(z^2/10)^2
- (x^2 + y^2 - z^2)
  *(-x^2 - 2*y^2 + 5*z^2)
"""
    ),

    (
        "5. Three ovals",
        """
(13*(x^2 - y^2))^2
- (35*x^2 + 35*y^2 - 4*z^2)
  *(3*x^2 + 3*y^2 + 16*z^2)
"""
    ),

    (
        "6. Four ovals",
        """
(z^2/10)^2
- (x^2 + 2*y^2 - z^2)
  *(-2*x^2 - y^2 + z^2)
"""
    ),
]


# ============================================================
#
# ============================================================

def automatic_input(prompt=""):
    print(prompt + "n")
    return "n"


old_input = builtins.input
builtins.input = automatic_input


# ============================================================
# ============================================================

try:

    for i, (name, quartic) in enumerate(quartics, start=1):

        print("\n")
        print("=" * 90)
        print(name)
        print("=" * 90)

        #
        quartic = " ".join(
            line.strip()
            for line in quartic.strip().splitlines()
        )

        # 
        replacement = "f=(" + quartic + ")"

        modified_source, number_replaced = re.subn(
            r'(?m)^\s*f\s*=.*$',
            replacement,
            vinzant_source,
            count=1
        )

        if number_replaced != 1:
            raise RuntimeError(
                "\n"
                "Check quartictype_fixed.sage 。"
            )

        # 
        temp_file = Path(
            "_vinzant_test_{:02d}.sage".format(i)
        )

        temp_file.write_text(modified_source)

        print("\nInput quartic:")
        print("f =", quartic)
        print("\nRunning Vinzant's original code...\n")

        try:
            load(str(temp_file))

        except SystemExit:
            # 
            pass

        finally:
            if temp_file.exists():
                temp_file.unlink()

finally:
    builtins.input = old_input


print("\n")
print("=" * 90)
print("Finished all six quartics.")
print("=" * 90)
