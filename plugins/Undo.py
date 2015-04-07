import pyrat


class Undo(pyrat.Worker):
    """
    Undo

    :author: Andreas Reigber
    """
    gui = {'menu': 'Tools', 'entry': 'Undo'}

    @classmethod
    def guirun(cls, viewer):
        if len(viewer.undo) >= 2:
            viewer.statusBar.setMessage(message=' Undo ', colour = 'R')

            actual = set(viewer.undo[-1][0])
            before = set(viewer.undo[-2][0])
            diff = actual.difference(before)

            pyrat.data.activateLayer(viewer.undo[-2][1])
            pyrat.data.delLayer(list(diff))

            viewer.undo.pop()
            viewer.undo.pop()

            viewer.statusBar.setMessage(message=' Ready ', colour='G')
            viewer.updateViewer()




