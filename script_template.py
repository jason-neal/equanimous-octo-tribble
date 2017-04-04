#!/usr/bin/python

# Basic argparse script template
import argparse


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Helpful description')
    parser.add_argument('pos_arg', help='A postitional argument', dtype=str)
    parser.add_argument('-f', "--flag", action="store_true", help='A toggle flag')
    args = parser.parse_args()
    return args


def main(pos_arg, flag=False):
    """Do main stuff."""
    pass


if __name__ == '__main__':
    args = vars(_parser())
    pos_arg = args.pop('pos_arg')  # positional arguments

    opts = {k: args[k] for k in args}  # keyword arguments

    main(pos_arg, **opts)
