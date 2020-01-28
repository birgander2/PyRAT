"""
AutoInstaller
=============

This script checks for missing dependencies/modules and tries to download them

"""
import platform
import sys
import importlib
from qgis.utils import iface


def checkInstallation(needed, condaPkgs, pipPkgs, qtModule="PyQt5.QtWidgets"):
    global _neededModules, _condaPackages, _pipPackages, _unknownModules
    global _qtModule
    _neededModules = needed
    _condaPackages = condaPkgs
    _pipPackages = pipPkgs
    _qtModule = None
    _unknownModules = []

    if not platform.python_version().startswith("3"):
        _print("ERROR: This python script requires Python 3.")
        sys.exit(1)

    for module in _neededModules:
        _testImport(module)

    if _unknownModules:
        _print("WARNING: Missing packages: " + str(_unknownModules))
        _testImport(qtModule)
        if qtModule not in _unknownModules:
            _qtModule = importlib.import_module(qtModule)
            _qtModule.QMessageBox.warning(
                    iface.mainWindow(),
                    "Missing packages",
                    "The following packages are missing: " +
                    str(_unknownModules))

        _testImport("conda")
        if "conda" in _unknownModules:
            # Anaconda not availible in current environment or not installed
            from subprocess import run, PIPE
            if(sys.platform == "linux" and
               run("conda", stdout=PIPE).returncode == 0):
                condaAvailible = 2
            elif sys.platform != "win32":
                try:
                    if(run("anaconda", stdout=PIPE).returncode == 0):
                        condaAvailible = 2
                    else:
                        condaAvailible = False
                except FileNotFoundError:
                    condaAvailible = False
            else:
                condaAvailible = False
        else:
            condaAvailible = 1

        _testImport("pip")
        pipAvailible = "pip" not in _unknownModules

        if condaPkgs and condaAvailible and _installWith("conda"):
            if condaAvailible == 1:
                _installCondaInternal()
            else:
                _installCondaSystem()
        elif pipPkgs and pipAvailible and _installWith(
                "pip (unstable, can corrupt your conda environment)"):
            _installPip()
        else:
            _print("Please install the missing packages "
                   "manually or create a new conda environment "
                   "with the requirements.yml file.")
            _print("Aborting.")
            if _qtModule:
                _qtModule.QMessageBox.critical(
                        iface.mainWindow(),
                        "Aborting",
                        "Please install the missing packages manually "
                        "or create a new conda environment "
                        "with the requirements.yml file.")
            sys.exit(1)

        _print("Installation of missing packages successful. Continuing.")
        if _qtModule:
            _qtModule.QMessageBox.information(
                    iface.mainWindow(),
                    "Installation successful",
                    "The installation of the missing packages was successful.")

    del _qtModule
    del _neededModules
    del _pipPackages
    del _condaPackages
    del _unknownModules


def _testImport(test):
    try:
        importlib.import_module(test)
    except ModuleNotFoundError:
        _unknownModules.append(test)


def _noWritePerms():
    _print("ERROR: The current user does not have "
           "write permissions to the target environment.")
    if _qtModule:
        _qtModule.QMessageBox.critical(
                iface.mainWindow(),
                "No write permissions",
                "You don't have write permissions to install the requested"
                " packages in the current python environment. Ask the "
                "environment administrator to install the packages or create "
                "a new conda environment with the requirements.yml file.")


def _unknownError():
    if _qtModule:
        _qtModule.QMessageBox.critical(
                iface.mainWindow(),
                "Aborting",
                "Errors occured during package installation.\r\n"
                "Please look into the console for more details.")


def _installCondaInternal():
    import conda
    import conda.cli.python_api as Conda
    try:
        _, _, return_code_int = Conda.run_command(
                Conda.Commands.INSTALL,
                *_condaPackages,
                stdout=sys.stdout,
                stderr=sys.stderr)
    except conda.CondaMultiError as e:
        if(len(e.errors) == 1 and
           type(e.errors[0]) is conda.exceptions.EnvironmentNotWritableError):
            _noWritePerms()
        else:
            _unknownError()
        raise e


def _installCondaSystem():
    from subprocess import run, PIPE
    p = None
    if(sys.platform == "linux"):
        p = run(["conda", "install", "-y"] + _condaPackages,
                universal_newlines=True, stderr=PIPE)
    else:
        # Windows
        p = run(["anaconda", "install", "-y"] + _condaPackages,
                universal_newlines=True, stderr=PIPE)
    if p.returncode != 0:
        _print(p.stderr)
        if("EnvironmentNotWritableError" in p.stderr or
           "Access is denied" in p.stderr):
            _noWritePerms()
        else:
            _unknownError()
        sys.exit(1)


def _print(string, *args, **kwargs):
    from qgis.core import QgsApplication
    QgsApplication.messageLog().logMessage(string)


def _installPip():
    class Cache():
        def __init__(self):
            self.log = str()

        def write(self, data):
            self.log += data

        def flush(self):
            pass
    cache = Cache()

    stderrcache = sys.stderr
    sys.stderr = cache
    if(sys.platform == "win32"):
        sysexeccache = sys.executable
        stdoutcache = sys.stdout
        sys.stdout = cache
        import os
        sys.executable = os.path.join(
                * os.path.split(sys.executable)[:-1], "python3.exe")

    import pip._internal
    code = pip._internal.main(
            ["install", "-q", "--upgrade-strategy",
             "only-if-needed"] + _pipPackages)

    sys.stderr = stderrcache
    if(sys.platform == "win32"):
        sys.stdout = stdoutcache
        sys.executable = sysexeccache

    if code != 0:
        _print(cache.log)
        if("[Errno 13] Permission denied" in cache.log or
           "Access is denied" in cache.log):
            _noWritePerms()
        else:
            _unknownError()
        sys.exit(1)


def _installWith(system):
    if _qtModule:
        ans = _qtModule.QMessageBox.question(
                iface.mainWindow(),
                "Install missing packages?",
                "Install missing packages with " + system + "?\r\n"
                "(You'll get a notification when finished)",
                _qtModule.QMessageBox.Yes | _qtModule.QMessageBox.No,
                _qtModule.QMessageBox.Yes)
        return ans == _qtModule.QMessageBox.Yes
    else:
        valid = {"": True, "yes": True, "y": True,
                 "no": False, "on": False, "n": False}
        while True:
            choice = input(
                    "Install missing packages with " + system + " [Y/n]? "
                    ).lower()
            if choice in valid:
                return valid[choice]
            else:
                _print("Please respond with 'yes' or 'no' (or 'y' or 'n').")
