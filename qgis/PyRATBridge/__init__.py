def classFactory(iface):
    """This function is needed for initalisation as QGIS-Plugin"""
    from .preLoader import PreLoader
    return PreLoader()
