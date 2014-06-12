def change_kwdefaults(func, **wrapkwargs):
    def inner(*args, **kwargs):
        for i in wrapkwargs.iteritems():  kwargs[i[0]] = i[1]
        return func(*args, **kwargs)
    return inner

    
    
    
    
    
    
    
    