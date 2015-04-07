import pyrat


class MyWizard(pyrat.WizardWorker):
    """
    Simple Wizard filter...
    """

    # gui = {'menu': 'Tools', 'entry': 'Template Wizard'}
    para = [
        {'var': 'win', 'value': 100, 'type': 'int',  'range': [3, 999],   'text': 'Window size', 'subtext': ['range', 'azimuth']},
        {'var': 'method', 'value': 'original', 'type': 'list', 'range': ['original', 'old', 'new'], 'text': 'Method'}
    ]

    def wizard(self, *args, **kwargs):
        pyrat.filter.complex2abs()
        pyrat.filter.boxcar(win=self.win)


def mywizard(*args, **kwargs):
    return MyWizard(*args, **kwargs).run(**kwargs)
