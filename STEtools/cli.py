import readline, sys, curses


def edit_input(prompt, default=''):
    """
    Let the user input a value, with an editable default value.
    
    :author: Internet
    """
    
    readline.set_startup_hook(lambda: readline.insert_text(default))
    try:
        return raw_input(prompt)
    finally:
        readline.set_startup_hook(None)

      
class ProgressBar():
    """
    Simple progress bar for the command line. 
   
    :author: Andreas Reigber
    """
    def __init__(self, message, max, **kwargs):
        curses.setupterm()
        terminal_width = curses.tigetnum('cols')
        if not sys.stdout.isatty(): terminal_width = 70
        self.message = message.ljust(25)
        if terminal_width < 50: self.message = self.message[:terminal_width/2]
        self.width   = terminal_width - len(self.message) - 10
        self.max     = max

    def __del__(self):
        print
        
    def update(self,val):
        """
        Updates the progress bar to the given value of progress
        """
        percent = float(val) / self.max
        hashes = '#' * int(round(percent * self.width))
        spaces = ' ' * (self.width - len(hashes))
        retline = "\r" if sys.stdout.isatty() else ""
        if sys.stdout.isatty() or val == 0:
            sys.stdout.write(retline+self.message+": [{0}] {1}%".format(hashes + spaces, int(round(percent * 100))))
            sys.stdout.flush()
           
