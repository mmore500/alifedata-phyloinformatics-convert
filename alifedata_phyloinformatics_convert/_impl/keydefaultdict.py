# adapted from https://github.com/Lyncs-API/lyncs.utils
# BSD 3-Clause License

# Copyright (c) 2020, Lyncs-API
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from collections import defaultdict
from contextlib import contextmanager


@contextmanager
def _setting(obj, attr, value, default=None):
    """
    Context manager that temporarily sets an attribute of an object.
    It yields the current value of the attribute if present,
    otherwise returns `default` (`None` by default).
    If the attribute did not exist, it is removed on exit,
    otherwise its original value is restored.
    """

    try:
        old = getattr(obj, attr)
        had = True
    except AttributeError:
        old = default
        had = False

    setattr(obj, attr, value)

    yield old

    if had:
        setattr(obj, attr, old)
    else:
        delattr(obj, attr)


class keydefaultdict(defaultdict):
    "A defaultdict that passes the key to the factory as argument"

    def __missing__(self, key):
        with _setting(self, "default_factory", lambda: fnc(key)) as fnc:
            return super().__missing__(key)
