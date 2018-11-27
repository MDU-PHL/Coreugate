#!/usr/bin/env python3

import logging
import argparse
import pathlib
from cleo import Application

from scripts.Coreugate import RunCoreugate

logging.basicConfig(filename='job.log',level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

application = Application(name = 'Coreugate', version= '2.0')
application.add(RunCoreugate())


def main():
    application.run()

if __name__ == '__main__':
    main()


