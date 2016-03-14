C++ coding rules
================

All ctools C++ code must be compatible with ISO standard C++98. Do not
use C++11 features (see `here <http://en.wikipedia.org/wiki/C%2B%2B#Standardization>`__
if you don't know what C++98 and C++11 is and note that C++03 is just
a minor correction of C++98).

The reason for this restriction is that C++ compilers / standard libraries
such as `gcc <http://gcc.gnu.org>`_ / `libstdc++ <http://gcc.gnu.org/libstdc++/>`_
or `clang <http://clang.llvm.org>`_ / `libc++ <http://libcxx.llvm.org>`_
only finished implementing C++11 in 2013 and using C++11 features in
ctools would make it impossible to use the code on many systems.


Code lay-out
^^^^^^^^^^^^

String quotes
^^^^^^^^^^^^^

Whitespace in Expressions and Statements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Comments
^^^^^^^^

Naming Conventions
^^^^^^^^^^^^^^^^^^                                

Programming Recommendations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

