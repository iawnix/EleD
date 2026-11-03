"""Public xTB adapter exceptions."""

class XTBError(RuntimeError):
    """Base class for xTB adapter errors."""


class XTBNativeNotFoundError(XTBError):
    """Raised when a shared libxtb cannot be located."""


class XTBRunError(XTBError):
    """Raised when xTB does not terminate normally."""



class DescriptorUnavailableError(ValueError):
    """A descriptor needs data not supplied by the selected calculation backend."""
