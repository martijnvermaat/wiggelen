"""
Helper functions for parsing wiggle tracks.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from collections import namedtuple


class ParseError(Exception):
    """
    Raised if a wiggle track cannot be parsed.
    """
    pass


# Line types.
class LineType(object):
    NONE, REGION, DATA = range(3)


# Track modes.
class Mode(object):
    VARIABLE, FIXED = range(2)


Data = namedtuple('Data', 'position span value')


# Create a new object encapsulating state for the parse function.
def create_state():
    return dict(mode=Mode.VARIABLE, span=1, start=None, step=None)


def parse(line, state):
    # Parse a line and return a tuple (line_type, data). The state dictionary is
    # modified and should be passed as such with the next call.
    if line.startswith('browser') or line.startswith('track') \
           or line.startswith('#') or line.isspace():
        # As far as I can see empty lines and comments are not allowed
        # by the spec, but I guess they exist in the real world.
        return LineType.NONE, None

    if line.startswith('variableStep'):
        try:
            fields = dict(map(lambda field: field.split('='),
                              line[len('variableStep'):].split()))
            state['mode'] = Mode.VARIABLE
            state['span'] = int(fields.get('span', 1))
            return LineType.REGION, fields['chrom']
        except (ValueError, KeyError):
            raise ParseError('Could not parse line: %s' % line)

    if line.startswith('fixedStep'):
        try:
            fields = dict(map(lambda field: field.split('='),
                              line[len('fixedStep'):].split()))
            state['mode'] = Mode.VARIABLE
            state['start'] = int(fields['start'])
            state['step'] = int(fields['step'])
            state['span'] = int(fields.get('span', 1))
            return LineType.REGION, fields['chrom']
        except (ValueError, KeyError):
            raise ParseError('Could not parse line: %s' % line)

    if state['mode'] == Mode.VARIABLE:
        try:
            position, value = line.split()
            return LineType.DATA, Data(
                position=int(position),
                span=state['span'],
                value=float(value) if '.' in value else int(value))
        except ValueError:
            raise ParseError('Could not parse line: %s' % line)

    if state['mode'] == Mode.FIXED:
        try:
            value = float(line) if '.' in line else int(line)
            state['start'] += state['step']
            return LineType.DATA, Data(
                position=start,
                span=state['span'],
                value=value)
        except ValueError:
            raise ParseError('Could not parse line: %s' % line)

    raise ParseError('Could not parse line: %s' % line)
