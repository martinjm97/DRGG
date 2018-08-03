from decimal import Decimal


def wrap_decimal(f):
    """ Converts all arguments to instances of Decimal for more precision """
    return lambda *args: f(*[Decimal(arg) for arg in args])
