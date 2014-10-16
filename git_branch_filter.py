#!/usr/bin/env python

"""
git filter to change branch=__BRANCH__ to reflect current branch.

This is useful for storing branch info in README.md, as well
as automatically merged travis.yml files, or other locations
where knowledge of the branch is useful.

"""

import re
import sys
import subprocess

SMUDGE = 'smudge'
CLEAN = 'clean'

branch = subprocess.check_output('git rev-parse --abbrev-ref HEAD'.split()).replace('\n', '')


def _error(msg):
    sys.stderr.write(msg+'\n')
    sys.stderr.flush()
    sys.exit(1)


def clean():
    for line in sys.stdin:
        if 'branch=__BRANCH__' in line:
            line = re.sub('branch=__BRANCH__', 'branch=%s' % branch, line)
        sys.stdout.write(line)


def smudge():
    for line in sys.stdin:
        if 'branch=%s' % branch in line:
            line = re.sub('branch=%s' % branch, 'branch=__BRANCH__', line)
        sys.stdout.write(line)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        _error('ERR: missing command line params %s or %s' % (SMUDGE, CLEAN))

    if sys.argv[1] == SMUDGE:
        smudge()
    elif sys.argv[1] == CLEAN:
        clean()
    else:
        _error('ERR: first arg must be %s or %s' % (SMUDGE, CLEAN))
