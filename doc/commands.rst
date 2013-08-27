Command line interface
======================

Some of the functionality in Wiggelen is provided through a simple command
line interface.

Since the average scientist is too lazy to write complete documentation,
you'll just find a quick dump of the command line help output below.

::

    martijn@hue:~$ wiggelen -h
    usage: wiggelen [-h]
                    {index,sort,scale,derivative,plot,coverage,merge,distance} ...

    Wiggelen command line interface.

    optional arguments:
      -h, --help            show this help message and exit

    subcommands:
      {index,sort,scale,derivative,plot,merge,distance}
                            subcommand help
        index               build index for wiggle track
        sort                sort wiggle track regions alphabetically
        scale               scale values in a wiggle track
        derivative          create derivative of a wiggle track
        plot                visualize wiggle tracks in a plot (requires
                            matplotlib)
        coverage            create coverage BED track of a wiggle track
        merge               merge any number of wiggle tracks in various ways
        distance            calculate the distance between wiggle tracks

Well, I guess nobody ever got fired for showing a quick example, so here you
go::

    martijn@hue:~$ wiggelen distance tests/data/*.wig
    A: tests/data/a.wig
    B: tests/data/b.wig
    C: tests/data/complex.wig
    D: tests/data/c.wig
    E: tests/data/empty.wig

          A     B     C     D     E
    A     x
    B   0.687   x
    C   0.000 0.687   x
    D   0.901 0.958 0.901   x
    E   0.974 0.952 0.974 0.748   x
