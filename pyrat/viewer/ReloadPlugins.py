import pyrat
import importlib


class ReloadPlugins(pyrat.Worker):
    """
    Reloads all plugins

    :author: Felix Weinmann
    """

    gui = {'menu': 'Tools', 'entry': 'Reload plugins'}

    def main(self, *args, **kwargs):
        importlib.reload(pyrat.plugins)

        if hasattr(pyrat, "app"):
            pyrat.app.menubar.clear()
            pyrat.app.makeMenu()
            pyrat.app.initPlugins()
            pyrat.app.makeActions()


@pyrat.docstringfrom(ReloadPlugins)
def reloadplugins(*args, **kwargs):
    return ReloadPlugins(*args, **kwargs).run(*args, **kwargs)
