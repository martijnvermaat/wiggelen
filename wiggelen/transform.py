"""
Various transformations on wiggle tracks.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

from collections import deque


# Difference directions.
_BACKWARD, _FORWARD, _CENTRAL = 0, 1, 2


def _divided_difference(walker, direction=_CENTRAL, step=None,
                        auto_step=False):
    queue_size = 3 if direction == _CENTRAL else 2
    current = 0 if direction == _FORWARD else 1
    queue = deque(maxlen=queue_size)

    for region, position, value in walker:
        # Oldest entry is automatically dropped when exceeding maxlen.
        queue.append((region, position, value))
        regions, positions, values = zip(*queue)

        # Initialize step size.
        if step is None and auto_step and len(queue) > 1:
            step = positions[1] - positions[0]

        # Make sure derivative is defined for current position.
        if (len(queue) != queue_size or
            regions[-1] != regions[0] or
            (step and any(b - a != step
                          for a, b in zip(positions, positions[1:])))):
            continue

        yield (regions[current], positions[current],
               (values[-1] - values[0]) / (positions[-1] - positions[0]))


def forward_divided_difference(walker, step=None, auto_step=False):
    """
    Derivative calculated by the forward divided difference method.

    .. note:: This transformation only works on walkers with numerical values.

    :arg walker: Generator yielding tuples of (region, position, value) per
        defined position.
    :type walker: generator(str, int, float)
    :arg step: Restrict calculation to positions that are this far apart
        (no restriction if `None`).
    :type step: int
    :arg auto_step: If `True` and `step=None`, automatically set `step` to a
        value based on the first two positions in `walker`.
    :type auto_step: bool

    :return: Tuple of (region, position, derivative value) per defined
        position in `walker` for which the derivative value is defined.
    :rtype: generator(str, int, float)
    """
    return _divided_difference(walker, direction=_FORWARD, step=step,
                               auto_step=auto_step)


def backward_divided_difference(walker, step=None, auto_step=False):
    """
    Derivative calculated by the backward divided difference method.

    .. note:: This transformation only works on walkers with numerical values.

    :arg walker: Generator yielding tuples of (region, position, value) per
        defined position.
    :type walker: generator(str, int, float)
    :arg step: Restrict calculation to positions that are this far apart
        (no restriction if `None`).
    :type step: int
    :arg auto_step: If `True` and `step=None`, automatically set `step` to a
        value based on the first two positions in `walker`.
    :type auto_step: bool

    :return: Tuple of (region, position, derivative value) per defined
        position in `walker` for which the derivative value is defined.
    :rtype: generator(str, int, float)
    """
    return _divided_difference(walker, direction=_BACKWARD, step=step,
                               auto_step=auto_step)


def central_divided_difference(walker, step=None):
    """
    Derivative calculated by the central divided difference method.

    .. note:: This transformation only works on walkers with numerical values.

    :arg walker: Generator yielding tuples of (region, position, value) per
        defined position.
    :type walker: generator(str, int, float)
    :arg step: Restrict calculation to positions that are this far apart. If
        `None`, automatically set `step` to a value based on the first two
        positions in `walker`.
    :type step: int

    :return: Tuple of (region, position, derivative value) per defined
        position in `walker` for which the derivative value is defined.
    :rtype: generator(str, int, float)
    """
    return _divided_difference(walker, direction=_CENTRAL, step=step,
                               auto_step=step is None)
