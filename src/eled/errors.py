"""Public xTB adapter exceptions."""

class XTBError(RuntimeError):
    """Base class for xTB adapter errors."""


class XTBNativeNotFoundError(XTBError):
    """Raised when a shared libxtb cannot be located."""


class XTBRunError(XTBError):
    """Raised when xTB does not terminate normally."""


class XTBEnvironmentError(XTBRunError):
    """Engine setup, allocation or extraction failed; do not skip as an SCF row."""


class XTBTimeoutError(XTBRunError):
    """An isolated calculation exceeded its wall-clock budget."""


class XTBWorkerError(XTBRunError):
    """An isolated calculation process exited without returning a result."""


class XTBCacheError(XTBError):
    """Persistent data is incompatible, corrupt or conflicts with this dataset."""



class DescriptorUnavailableError(ValueError):
    """A descriptor needs data not supplied by the selected calculation backend."""
