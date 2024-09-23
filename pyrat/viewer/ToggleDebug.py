import pyrat


class ToggleDebug(pyrat.Worker):
    """
    This class toggles the debug mode behaviour.

    :author: Felix Weinmann
    """

    gui = {'menu': 'Tools', 'entry': 'Toggle debug mode'}

    def main(self, *args, **kwargs):
        if pyrat._debug:
            pyrat._debug = False
            pyrat.pyrat_debug(False)
        else:
            pyrat._debug = True
            pyrat.pyrat_debug(True)


@pyrat.docstringfrom(ToggleDebug)
def toggledebug(*args, **kwargs):
    return ToggleDebug(*args, **kwargs).run(*args, **kwargs)
