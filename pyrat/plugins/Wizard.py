import pyrat


class TestGroup(pyrat.GroupWorker):
    """
    Example for a simple processing group...
    """

    # gui = {'menu': 'Tools', 'entry': 'Template Wizard'}
    para = [
        {'var': 'win', 'value': [10, 10], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
         'subtext': ['range', 'azimuth']},
        {'var': 'method', 'value': 'original', 'type': 'list', 'range': ['original', 'old', 'new'], 'text': 'Method'}
    ]

    def group(self, *args, **kwargs):
        pyrat.filter.complex2abs()
        pyrat.filter.boxcar(win=self.win)


@pyrat.docstringfrom(TestGroup)
def testgroup(*args, **kwargs):
    return TestGroup(*args, **kwargs).run(*args, **kwargs)
