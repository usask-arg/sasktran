from sasktranif._sasktranif import functionfail
import functools


class SasktranError(functionfail):
    pass


def wrap_skif_functionfail(func, message=None):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        skif_exception = None
        try:
            return func(*args, **kwargs)
        except functionfail as e:
            # Store the skif exception, but don't raise here or else we will get funny chained exceptions
            skif_exception = e

        raise SasktranError(skif_exception)
    return wrapper
