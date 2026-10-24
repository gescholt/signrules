"""
Command line of the extension screen (``signrules.screen``): the candidates for minimal impossible
n-species patterns by the deletion lemma, and the sign rules on each.

    python scripts/extension_screen.py 5 --workers 12                  # from the n=4 census: 83 candidates
    python scripts/extension_screen.py 6 --base known --workers 12     # from the 24 n=5 kills: 143 candidates
    python scripts/extension_screen.py 6 --candidates-only             # count candidates only

Equivalent to the console script ``signrules-screen``.
"""
from signrules.screen import main

if __name__ == "__main__":
    main()
