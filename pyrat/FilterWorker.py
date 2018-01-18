import pyrat
import logging


class FilterWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(FilterWorker, self).__init__(*args, **kwargs)

    def main(self, *args, **kwargs):
        """
        Main routine doing the (parallel) block processing, calling the (overloaded) filter method.
        Before the actual processing, pre() is called, as well as post() after completing it.
        """

        if self.checkinput():
            self.pre(*args, **kwargs)
            newlayer = self.layer_process(self.filter, silent=False, **kwargs)
            pyrat.data.activateLayer(newlayer)
            self.post(*args, **kwargs)

            if self.delete is True:
                pyrat.delete(self.input)
            return newlayer
        else:
            return False

    def filter(self, data, *args, **kwargs):
        """
        The actual filter routine (to be overloaded)
        """
        logging.error(self.name + ': No filter method defined')
        return False

    def pre(self, *args, **kwargs):
        """
        The preprocessing routine (to be overloaded)
        """
        pass

    def post(self, *args, **kwargs):
        """
        The postprocessing routine (to be overloaded)
        """
        pass

    def info(self):
        print(self.__doc__)


