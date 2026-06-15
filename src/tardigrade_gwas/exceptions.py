"""Package-specific exception hierarchy."""

class TardigradeError(Exception):
    """Base exception for Tardigrade-specific failures."""

class InputTableError(TardigradeError, ValueError):
    """Raised when an input table cannot be interpreted."""
