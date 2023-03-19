import string
import typing


def parse_ancestor_list(raw: str) -> typing.List[int]:
    assert not any(
        c in raw for c in string.whitespace
    ), "Whitespace separated ancestor list not supported."
    try:
        ancestor_list = eval(raw)
        if ancestor_list == [None]:
            return []
        else:
            assert None not in ancestor_list
            return ancestor_list
    except NameError:
        # if ancestor list contains a placeholder string like NONE,
        # it will trigger a NameError when eval'ed
        return []
