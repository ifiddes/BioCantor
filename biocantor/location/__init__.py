"""
:class:`Location` objects represent features defined with respect to a coordinate system. The :class:`Location` API
includes rich feature arithmetic methods. Additionally, a :class:`Location` object can define its relationship with
a :class:`Parent` object, situating it within a potentially arbitrary hierarchy of coordinate systems; the
:class:`Location` API provides rich coordinate and location conversion methods.
"""

from biocantor.location.location import Location
from biocantor.location.strand import Strand
from biocantor.parent import _make_parent_dispatch, Parent
from biocantor.location.location_impl import SingleInterval, CompoundInterval, EmptyLocation  # noqa F401


@_make_parent_dispatch.register(Location)
def _(obj) -> Parent:
    return Parent(location=obj)


@_make_parent_dispatch.register(Strand)
def _(obj) -> Parent:
    return Parent(strand=obj)
