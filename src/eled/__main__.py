"""Allow ``python -m eled doctor`` in the selected Python environment."""

from .diagnostics import main

raise SystemExit(main())
