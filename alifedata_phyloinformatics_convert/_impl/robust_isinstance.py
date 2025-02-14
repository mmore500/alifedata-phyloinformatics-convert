def robust_isinstance(obj: object, classinfo: object) -> bool:
    try:
        return isinstance(obj, classinfo)
    except TypeError:
        return False