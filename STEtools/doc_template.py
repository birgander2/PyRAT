"""A template for Python docstrings for Sphinx (sphinx-apidoc).

Here are a few tentative conventions that this example uses:

- Docstrings in functions are self-explanitory and should always
  indicate:

  * **name** (via :param...:) and **type** (via :type...:)of each
      parameter and the **author** (:author:...)
  * the **return** value

- The docstring after the class declaration should indicate:

  * author(s)
  * public **member variables** that might be of interest to others.
      These are declared using :ivar...: and the description should
      indicate the type of the member variable.
  * **name** and **type** of parameters to the constructor (as above).

- The docstring of class member functions does not need to include the author.

- When referencing classes that we have definied ourselves
  (e.g. as a parameter **type**), include a hyperlink reference. For example,
  the docstring of ClassB in this module references ClassA using 
  :py:class:`~Common.doc_template.ClassA`.

Note that all docstrings for Sphinx use so-called reStructuredText
(e.g. you can use **<something>** to make <something> bold in the online
documentation).

"""


def new_funct(arg1, opt1=True):
    """ Increments arg1 if opt1 is False.

    :author: Marc Jaeger
    
    :param arg1: The integer to increment
    :type arg1: int
    :param opt1: Increment flag
    :type opt1: bool

    :returns: arg1 incremented.

    """
    if (opt1 != True):
        return arg1+1
    return arg1



class ClassA(object):
    """Provides incrementing functionality

    :author: Marc Jaeger

    :ivar count: A counter.
    :param count: Initial counter value.
    :type count: int

    """

    def __init__(self,count):
        self.count = count


    def increment(self):
        """Increments the counter by one.

        :returns: The new value of the counter.

        """
        self.count += 1
        return self.count




class ClassB(ClassA):
    """Extends :py:class:`~Common.doc_template.ClassA`
    to support multiple incrementation
    
    :author: Marc Jaeger    

    :param count: Initial counter value.
    :type count: int

    """

    def __init__(self,count):
        super(ClassB,self).__init__(count)
        
    def incrementN(self,n):
        """Increments the class counter n times

        :param n: The number of times to increment
        :type n: int

        """
        for m in range(n):
            self.increment()

        return self.count



def testTemplate():
    """Tests template functions/classes

    :author: Marc Jaeger
    
    """

    print 'new_func: ({},{})'.format(new_funct(42),new_funct(42,False))

    a = ClassA(41)
    print 'ClassA: ({},{})'.format(a.increment(),a.increment())

    b = ClassB(35)
    print 'ClassB: {}'.format(b.incrementN(7))



if __name__ == '__main__':
    testTemplate()
