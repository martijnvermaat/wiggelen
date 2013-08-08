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

    # As an optimization, we first check for the common case of a line with
    # data. It must always start with a number (either position or value).
    if line[0] in '0123456789.':
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
                position = state['start']
                state['start'] += state['step']
                return LineType.DATA, Data(
                    position=position,
                    span=state['span'],
                    value=float(line) if '.' in line else int(line))
            except ValueError:
                raise ParseError('Could not parse line: %s' % line)

    if line[:7] == 'browser' or line[:5] == 'track' \
           or line[0] == '#' or line in ('\n', '\r\n', '\r'):
        # As far as I can see empty lines and comments are not allowed
        # by the spec, but I guess they exist in the real world.
        return LineType.NONE, None

    if line[:12] == 'variableStep':
        try:
            fields = dict(map(lambda field: field.split('='),
                              line[len('variableStep'):].split()))
            state['mode'] = Mode.VARIABLE
            state['span'] = int(fields.get('span', 1))
            return LineType.REGION, fields['chrom']
        except (ValueError, KeyError):
            raise ParseError('Could not parse line: %s' % line)

    if line[:9] == 'fixedStep':
        try:
            fields = dict(map(lambda field: field.split('='),
                              line[len('fixedStep'):].split()))
            state['mode'] = Mode.FIXED
            state['start'] = int(fields['start'])
            state['span'] = int(fields.get('span', 1))
            # Though not valid by the specification, we accept fixedStep
            # definitions without a step argument (it's also accepted by the
            # UCSC Genome Browser).
            # Issue: https://github.com/martijnvermaat/wiggelen/issues/1
            state['step'] = int(fields.get('step', min(1, state['span'])))
            return LineType.REGION, fields['chrom']
        except (ValueError, KeyError):
            raise ParseError('Could not parse line: %s' % line)

    raise ParseError('Could not parse line: %s' % line)
