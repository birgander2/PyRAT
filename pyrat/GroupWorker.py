import pyrat
import logging


class GroupWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(GroupWorker, self).__init__(*args, **kwargs)

    def main(self, *args, **kwargs):
        if self.checkinput():
            self.group(*args, **kwargs)
            newlayer = pyrat.data.active
            if len(newlayer) == 1:
                newlayer = newlayer[0]
            return newlayer
        else:
            return False

    def group(self, *args, **kwargs):
        print("ERROR: group() method not overloaded")
        return False, False
