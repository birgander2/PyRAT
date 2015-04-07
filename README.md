# j-branch of PyRat (Python RadarTools)

## Installation and git usage

Create a directory:

    mkdir ~/progs/pyrat
    cd ~/progs/pyrat

Clone git reposotory:

    git clone https://mxn@bitbucket.org/mxn/pyrat.git

Update (pull) repository:

    git pull

Current status:

    git status

Adding new files to repository (only if this file is not yet tracked by git):

    git add filename

Commit & push changes:

    git commit –a –m "Useful doc-string describing the change"
    git push

In order to be able use the PyRat module, update your PYTHONPATH (write e.g. in .bashrc):

    export PYTHONPATH=~/progs/pyrat:$PYTHONPATH

## Content


* PyRat/ folder contains the main PyRat module, including main classes,
  library functions, and tests

* docs/ folder contains documentation

* icons/ folder contains the artwork for various icons of the GUI application
