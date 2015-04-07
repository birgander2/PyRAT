== COMPILE / BUILD ==
python setup.py build_ext --inplace

== INSTALL ==
# with root rights
python setup.py install

# otherwise / for user only
python setup.py install --user

== PACKAGE ==
python setup.py sdist

