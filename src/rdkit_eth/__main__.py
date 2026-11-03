"""Allow ``python -m rdkit_eth doctor`` in the selected Python environment."""

from .diagnostics import main

raise SystemExit(main())
