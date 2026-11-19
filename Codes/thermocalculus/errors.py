"""Exception types raised by the engine. Every failure mode has its own class
so that scenarios can report *why* a quantity could not be computed."""


class ThermoError(Exception):
    """Base class for all engine errors."""


class UnitError(ThermoError):
    """Unknown unit, missing unit, or a unit of the wrong dimension."""


class MissingDataError(ThermoError):
    """A target cannot be computed from the data supplied."""

    def __init__(self, message, missing=()):
        super().__init__(message)
        self.missing = list(missing)


class InconsistentInputError(ThermoError):
    """Supplied values contradict each other or violate a stability condition."""


class OutOfRangeError(ThermoError):
    """A table or correlation was evaluated outside its range of validity."""


class ConvergenceError(ThermoError):
    """An iterative solve failed."""


class InputError(ThermoError):
    """A scenario file is malformed."""


class SingularDerivativeError(ThermoError):
    """(dX/dY)_Z does not exist at this state because Y cannot be varied at
    constant Z (for example (dP/dT)_H in an ideal gas, where H depends on T only)."""
