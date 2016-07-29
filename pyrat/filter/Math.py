from __future__ import print_function
import numexpr as ne
import pyrat


class MathExpr(pyrat.FilterWorker):
    """
    Arithmetic layer caluculator analoguous (and using) numexpr. The first layer provided
    is assigned to A, the second to B, and so on. Result is calculated according to
    a mathematical string expression

    Example:
    filter.mathexpr(expr="A+B*B", layer=[lay1, lay2])
    """
    para = [
        {'var': 'expr', 'value': '', 'type': 'str', 'text': 'Expression'}
    ]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):

        letter = 'A'
        for var in array:
            vars()[letter] = var
            letter = chr(ord(letter) + 1)
        out = ne.evaluate(self.expr)
        return out


@pyrat.docstringfrom(MathExpr)
def mathexpr(*args, **kwargs):
    return MathExpr(*args, **kwargs).run(**kwargs)
