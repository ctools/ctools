Python coding rules
===================

The Python coding rules have been inspired from the Style Guide for
Python Code `PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_.
You can check whether you code complies to these rules using the
`pep8 tool <https://github.com/jcrocholl/pep8>`_ that you should run
on your Python code before committing.


Dependencies
^^^^^^^^^^^^

- All ctools code must be compatible with Python 2.3 or higher, and
  Python 3.0 and higher.

- ctools should be importable with no dependencies other than the
  `Python Standard Library <https://docs.python.org/release/2.3/lib/lib.html>`_.


Code lay-out
^^^^^^^^^^^^

- Put imports on separate lines, e.g.:

  .. code-block:: python

     Yes: import os
          import sys

     No:  import sys, os

- Put imports at the top of the file, just after any module comments
  and docstrings, and before module globals and constants.

- Use absolute imports, e.g.:

  .. code-block:: python

     import mypkg.sibling
     from mypkg import sibling
     from mypkg.sibling import example


String quotes
^^^^^^^^^^^^^

- Use single quotes for strings so that the string can contain double
  quote characters which avoid backslashes, e.g.:

  .. code-block:: python
 
     Yes: print('The Parameter "infile" is missing')

     No:  print("The Parameter \"infile\" is missing")


Whitespace in Expressions and Statements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Do not put extraneous whitespace immediately inside parentheses,
  brackets or braces, e.g.:

  .. code-block:: python
 
     Yes: spam(ham[1], {eggs: 2})

     No:  spam( ham[ 1 ], { eggs: 2 } )

- Do not put extraneous whitespace immediately before a comma, semicolon,
  or colon, e.g.:

  .. code-block:: python
 
     Yes: if x == 4: print x, y; x, y = y, x

     No:  if x == 4 : print x , y ; x , y = y , x

- Do not put extraneous whitespace immediately before the open parenthesis
  that starts the argument list of a function call, e.g.:

  .. code-block:: python
 
     Yes: spam(1)

     No:  spam (1)

- Do not put extraneous whitespace immediately before the open parenthesis
  that starts an indexing or slicing, e.g.:

  .. code-block:: python
 
     Yes: dct['key'] = lst[index]

     No:  dct ['key'] = lst [index]


Comments
^^^^^^^^

- Always make a priority of keeping the comments up-to-date when the code
  changes! Comments that contradict the code are worse than no comments. 

- Comments should be complete sentences. If a comment is a phrase or
  sentence, its first word should be capitalized, unless it is an identifier
  that begins with a lower case letter (never alter the case of identifiers!).

- If a comment is short, the period at the end can be omitted. Block
  comments generally consist of one or more paragraphs built out of complete
  sentences, and each sentence should end in a period.

- Follow `Strunk & White <https://en.wikipedia.org/wiki/The_Elements_of_Style>`_.
  Make every word tell. Omit needless words. Use active voice. Use parallel
  construction on parallel concepts.

- Block comments generally apply to some (or all) code that follows them,
  and are indented to the same level as that code. Each line of a block
  comment starts with a # and a single space (unless it is indented text
  inside the comment). Paragraphs inside a block comment are separated by
  a line containing a single #, e.g.:

  .. code-block:: python
 
     # Compute the energy boundaries.
     #
     # The energy boundaries are computed from the lower and upper energy
     # thresholds that are stored in the effective area components of the
     # response functions.
     ebounds = gammalib.GEbounds()
     ...

- Use inline comments sparingly. Use meaningful variable names instead.
  But sometimes, inline comments can be useful, e.g.:

  .. code-block:: python

     x = x + 1  # Compensate for border


User Documentation
^^^^^^^^^^^^^^^^^^

- Write documentation strings (a.k.a. docstrings) for all modules, functions,
  classes, and public methods.
  Read `PEP 257 <http://www.python.org/dev/peps/pep-0257/>`_ to learn
  about the general Python conventions for writing good documentation
  strings.

- The docstring conventions follow the common standards developed by
  `NumPy <http://www.numpy.org>`_ and `SciPy <http://www.scipy.org>`_
  and that are also adopted by
  `AstroPy <http://docs.astropy.org/en/stable/development/docrules.html>`_.
  This docstring standard uses
  `re-structured text (reST) <http://docutils.sourceforge.net/rst.html>`_
  syntax and is rendered using
  `Sphinx <http://www.sphinx-doc.org/en/stable/>`_.
  While a rich set of markup is available, limit yourself to a very basic
  subset, in order to provide docstrings that are easy to read on text-only
  terminals.

- Format each docstring as follows (sections can be omitted if they do
  not apply, for example you do not need to specify that nothing is 
  returned):

  .. code-block:: python

     """
     Extract calibrations for a mission from a calibration database.

     .. note:: Deprecated in ctools 1.1
               The `instrument` string will be removed in ctools 2.0 because
               it's in fact not used by the method.

     The method extract calibration names from a calibration database by
     scanning all index files.

     Parameters
     ----------
     caldb :  `~gammalib.GCaldb`
         Calibration database.
     instrument : str
         Instrument names.
     debug : bool, optional
         Turn on debugging if needed.

     Returns
     -------
     missions : list of str
         List of mission names.

     Raises
     ------
     ValueError
         If `instrument` is not recognised as among the ones implemented.

     See Also
     --------
     irfs : Return all irfs for a given calibration.

     Notes
     -----
     See https//my.site.edu for documentation. The method computes

     .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

     where the value of :math:`\omega` is larger than 5.

     References
     ----------
     .. [1] J. Public, "The art of programming", Computers & Geosciences,
        vol. 22, pp. 585-588, 1996.

     Examples
     --------
     >>> calibrations(caldb, "CTA", debug=False)
     ['Prod2', 'Prod3']
     """

- For the parameter types, be as precise as possible. Below are a few examples
  of parameters and their types.

  .. code-block:: python

     Parameters
     ----------
     filename : str
     copy : bool
     dtype : data-type
     iterable : iterable object
     shape : int or tuple of int
     files : list of str

- For classes, use the same sections as outlined above. The constructor
  (``__init__``) should also be documented here, the Parameters section of
  the docstring details the constructors parameters. An Attributes section,
  located below the Parameters section, may be used to describe class
  variables:

  .. code-block:: python

     Attributes
     ----------
     x : float
         The X coordinate.
     y : float
         The Y coordinate.

- In general, it is not necessary to list class methods. In some cases,
  however, a class may have many methods, and then it becomes useful to have
  an additional ``Methods`` section:

  .. code-block:: python

    Methods
    -------
    read(filename)
        Read a table from a file
    sort(column, order='ascending')
        Sort by `column`


Naming Conventions
^^^^^^^^^^^^^^^^^^

- Use short lowercase abbreviate words for cscripts, e.g.:

  .. code-block:: python

     cspull
     csspec
     cslightcrv

- Use ``lower_case_with_underscores`` for functions, methods, and variables.

- Use ``self`` for the first argument to instance methods, e.g.:

  .. code-block:: python

     def __init__(self, name):
        self._name = name

- Use one leading underscore for private methods and instance variables,
  e.g.:

  .. code-block:: python

     self._get_energy_boundaries()
     self._has_ebounds = True

- Always decide whether a class's methods and instance variables
  (collectively: "attributes") should be public or private. If in doubt,
  choose private; it's easier to make it public later than to make
  a public attribute private.


.. To be written:

   Programming Recommendations
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^
