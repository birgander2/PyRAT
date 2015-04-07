import pyrat
import code


class Console1(pyrat.Worker):
    """
    Python Console to be called from GUI
    """

    gui = {'menu': 'Tools', 'entry': 'Python console'}

    @classmethod
    def guirun(cls, array, *args, **kwargs):
        variables = globals().copy()
        variables.update(locals())

        shell = code.InteractiveConsole(variables)
        shell.push("from pyrat import *")
        shell.interact(" \n <Ctrl>-D to exit\n")



