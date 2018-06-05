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


Dependencies
^^^^^^^^^^^^

- All ctools code must only depend on the ISO standard C++98 libraries and the
  GammaLib library.


Code lay-out
^^^^^^^^^^^^

- Each function or method starts with a curly bracket in the line following
  the function or method name. The return type is in the same line as the
  function or method name, e.g.:

  .. code-block:: cpp

     int function(void)
     {
         int i = 0;
         ...
         return i;
     }

- Code blocks should be encompassed in curly brackets, even if the block
  consists only of a single line.

- The opening curly bracket of a block starts in the same line as the
  related statement, e.g.:

  .. code-block:: cpp

      for (int i = 0; i < 10; ++i) {
          sum += i;
      }

- Separate code elements by spaces, e.g.:

  .. code-block:: cpp

     Yes: int i = 0;

     No:  int i=0;

- Align successive class definition on the member function name, e.g.:

  .. code-block:: cpp

     void        log10GeV(const double& eng);
     void        log10TeV(const double& eng);
     std::string print(void) const;


.. To be written:

   String quotes
   ^^^^^^^^^^^^^

   Whitespace in Expressions and Statements
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   Comments
   ^^^^^^^^


Naming Conventions
^^^^^^^^^^^^^^^^^^

- Prefix data members with ``m_``, e.g.:

  .. code-block:: cpp

     m_num
     m_response
     m_grid_length
     m_axis_dir_qual


Programming Recommendations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Each class shall have the protected methods ``init_members``,
  ``copy_members``, and ``free_members`` to handle initialisation,
  copying, and eventually freeing of allocated memory for data members:

  .. code-block:: cpp

     // Initialise data members
     void Class::init_members(void)
     {
         m_elements = 10;
         m_array    = new double[m_elements];
         for (int i = 0; i < m_elements; ++i) {
             m_array[i] = 0.0;
         }
         ...
     }

     // Copy data members
     void Class::copy_members(const &Class class)
     {
         m_elements = class.m_elements;
         m_array    = new double[m_elements];
         for (int i = 0; i < m_elements; ++i) {
             m_array[i] = class.m_array[i];
         }
         ...
     }

     // Free data members
     void Class::free_members(void)
     {
         if (m_array != NULL) delete [] m_array;
         ...
     }

- Each class shall have at least a void constructor, a copy constructor,
  a destructor and an assignment operator. Additional constructors and
  operators can be implemented as required. The following example shows
  the basic implementation for these 4 methods. Due to the usage of the
  ``init_members``, ``copy_members``, and ``free_members``, most classes
  will have exactly this kind of syntax:

  .. code-block:: cpp

     // Void constructor
     ctnice::ctnice(void)
     {
         init_members();
         return;
     }

     // Copy constructor
     ctnice::ctnice(const ctnice& nice)
     {
         init_members();
         copy_members(nice);
         return;
     }

     // Destructor
     ctnice::~ctnice(void)
     {
         free_members();
         return;
     }

     // Assignment operator
     ctnice& ctnice::operator=(const ctnice& nice)
     {
         if (this != &nice) {
             free_members();
             init_members();
             copy_members(nice);
         }
         return *this;
     }

- Do not use macros.

- Do not use ``#define`` directives for the declaration of constants. Use
  ``const`` instead.

- Do not use ``std::strncpy``, ``std::memcpy`` or similar as these functions
  are corrupted on some systems.

- If possible, pass arguments by reference.

- Output arguments should be passed as pointers.

- Use C++ (``std::string``) instead of C-style (``char*``) strings.

- Use C++ casts instead of C-style casts.

- Avoid using templates.

- Do not use an integer for a floating point argument (i.e. write 10.0
  instead of 10). Some older compilers give an error when using
  integers in some floating point functions, such as log10().

- Where possible (and appropriate), use ``std::vector`` containers instead
  of allocating memory. In other words: avoid direct memory allocation with
  ``new``.

- Use the ``std::`` namespace prefix where possible; write for example

  .. code-block:: cpp

     std::sin(angle);
     std::cos(angle);

  You may not believe it, but droping the ``std::`` may on some systems
  lead to considerably slower code for trigonometric functions!

- Use ``explicit`` for constructors with single arguments to prevent
  unintended type conversions. The only exception to this rule is the
  copy constructor or type conversion constructors.

- Specify ``void`` for function or method definitions without arguments,
  e.g.:

  .. code-block:: cpp

     Yes: void function(void)

     No:  void function()

- Use pre-incrementation in loops (pre-incrementation is faster than
  post-incrementation), e.g.:

  .. code-block:: cpp

     for (int i = 0; i < 10; ++i) {
         sum += i;
     }
