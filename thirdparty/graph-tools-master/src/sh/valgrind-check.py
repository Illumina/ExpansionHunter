#!/usr/bin/env python
# Simple checker for whether valgrind found errors

import sys
import xml.etree.ElementTree as ElementTree

e = ElementTree.parse(sys.argv[1])

states = [x.find('state').text for x in e.findall('status')]
errors = [x.find('kind').text for x in e.findall('error')]

if "RUNNING" not in states or "FINISHED" not in states:
    raise Exception("Valgrind didn't run successfully, states seen: %s" % str(states))

if errors:
    raise Exception("Valgrind found some errors in %s: %s" % (sys.argv[1], str(errors)))

sys.exit(0)
