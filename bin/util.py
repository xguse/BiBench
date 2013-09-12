###################################################################################################
###---------------------------------------------------------------------------------------------###
### This file is part of bibench (Biclustering Benchmarking)			                        ###
### Copyright (c) 2011,                                                                         ###
### By:    Kemal Eren,                                                                          ###
##         Mehmet Deveci,                                                                       ###
###        Onur Kucuktunc,                                                                      ###
###        Umit V. Catalyurek                                                                   ###
###                                                                                             ###
###---------------------------------------------------------------------------------------------###
### For license info, please see the README.txt and LICENSE.txt files in the main directory.    ###
###---------------------------------------------------------------------------------------------###
###################################################################################################
"""
Utility functions for generate_dataset.py and run_algorithm.py.

"""

import inspect

def get_nargs(arg, mult):
    nargs = None
    if not mult is None:
        if arg in mult:
            nargs = '+'
    return nargs


def extract_help(string):
    """
    Given the docstring of a function, tries to extract just the
    summary.

    Assumes the summary is in the first lines before a double newline.

    """

    lines = [l.strip() for l in string.split('\n')]
    try:
        while lines[0] == '':
            lines.pop(0)
    except IndexError:
        return ''
    helpstring = []
    if lines.count('') == 0:
        helpstring = lines
    else:
        while lines[0] != '':
            helpstring.append(lines.pop(0))
    return '\n'.join(helpstring)


#automatically create parsers for each type of dataset
def create_subparsers(add_parser, cmds, parents, ignore=None, mult=None):
    """
    Call add a subparser for each command, grabbing its arguments
    from the associated function definition.

    Args:
        * add_parser: The add_parser() function of a argparse._SubParsersAction
        * cmds: A dictionary in which keys are strings and values are functions.
            The keys are turned into subparser commands, and the associated function
            is used to populate its arguments.
        * parents: A list of ArgumentParser instances to be parents for each
            subparser
        * mult: If provided, a list of arguments for which nargs='+'.

    """
    for cmd, f in cmds.iteritems():
        help = extract_help(f.func_doc)
        subparser = add_parser(cmd,
                               help=help,
                               parents=parents)
        subparser.set_defaults(func=f)

        argspec = inspect.getargspec(f)
        firstdefault = -len(argspec.defaults)
        required = argspec.args[:firstdefault]
        optional = zip(argspec.args[firstdefault:], argspec.defaults)
        for r in required:
            if r in ignore:
                continue
            subparser.add_argument(r, nargs=get_nargs(r, mult))
        for o, d in optional:
            if o in ignore:
                continue
            subparser.add_argument('--' + o,
                                   default=d,
                                   nargs=get_nargs(o,mult))


def guess_type(arg):
    """
    Try to convert an arg to a number, or list of numbers.

    Warning: For lists with only one member, returns that member,
    converted.

    """
    if getattr(arg, '__iter__', False):
        if len(arg) == 1: #assume that single arg should not be a list.
            arg = arg[0]
        else:
            return [guess_type(a) for a in arg]
    try:
        arg = int(arg)
    except:
        try:
            arg = float(arg)
        except (ValueError, TypeError):
            pass
    return arg


def guess_args_type(kwargs):
    return {k : guess_type(v) for k, v in kwargs.iteritems()}


def get_kwargs(parser):
    import sys
    kwargs = vars(parser.parse_args())
    kwargs = {k:v for k, v in kwargs.iteritems() if not v is None}
    return guess_args_type(kwargs)
