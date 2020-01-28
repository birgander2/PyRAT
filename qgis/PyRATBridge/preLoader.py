import os
from qgis.utils import iface
from qgis.PyQt.QtWidgets import QMessageBox, QFileDialog


class PreLoader:

    def unload(self):
        if hasattr(self, "mainPlugin"):
            self.mainPlugin.unload()

    def initGui(self):
        from .autoInstaller import checkInstallation
        checkInstallation(
                ["numexpr", "scipy", "h5py", "PIL", "Cython", "matplotlib",
                 "skimage", "mako", "pyqtgraph", "pywt", "mkl",
                 "lxml"],
                ["-c", "defaults", "-c", "conda-forge", "numexpr", "scipy",
                 "h5py", "pillow", "cython", "matplotlib", "scikit-image",
                 "mako", "pyqtgraph>=0.10", "cvxpy", "pywavelets",
                 "mkl-service", "hdf5=1.10.4", "lxml"],
                ["numexpr", "scipy", "h5py", "Pillow", "Cython", "matplotlib",
                 "scikit-image", "Mako==1.0.12", "pyqtgraph", "PyWavelets",
                 "mkl", "lxml"],
                qtModule="qgis.PyQt.QtWidgets")

        import sys
        crashFile = os.path.join(os.path.dirname(__file__), "CHECK_CRASH")
        if(os.path.exists(crashFile)):
            if(sys.platform == "linux" and
               "CONDA_PREFIX" in os.environ and
               "LD_PRELOAD" not in os.environ):
                ans = QMessageBox.question(
                    iface.mainWindow(),
                    "Crash on last startup",
                    "The loading of PyRAT seems to have produced a crash.\n\r"
                    "The issue might be solved with a "
                    "LD_PRELOAD=$CONDA_PATH/lib/libedit.so statement. "
                    "Do you want to add this now?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.Yes)
                if(ans == QMessageBox.Yes):
                    _writeLDPRELOAD()
                    os.remove(crashFile)
                    return
            else:
                QMessageBox.information(
                            iface.mainWindow(),
                            "Crash on last startup",
                            "PyRAT seems to have produced a crash on last "
                            "startup. Please report the issue if the crash "
                            "continues to occur.")

        self.writeFail = False
        try:
            with open(crashFile, "w") as file:
                file.write("Checking pyrat loading...")
        except Exception:
            self.writeFail = True

        try:
            import pyrat
        except ModuleNotFoundError:
            print("WARNING: No PyRAT installation found")
            QMessageBox.warning(iface.mainWindow(),
                                "PyRAT not found",
                                "PyRAT was not found in the current "
                                "python environment.")
            while(True):
                path = QFileDialog.getExistingDirectory(
                        iface.mainWindow(),
                        "PyRAT-Folder")
                if path == "":
                    return
                if not os.path.isdir(os.path.join(path, "pyrat")):
                    QMessageBox.warning(
                            iface.mainWindow(),
                            "No PyRAT installation",
                            "The selected folder contains no pyrat-Module."
                            " Please select a valid PyRAT folder. "
                            "If you currently don't have PyRAT: "
                            "https://github.com/birgander2/PyRAT")
                else:
                    import shutil
                    shutil.copytree(
                            os.path.join(path, "pyrat"),
                            os.path.join(sys.path[0], "pyrat"))
                    break

        if not self.writeFail:
            os.remove(crashFile)

        from .mainPlugin import PyRATBridge
        self.mainPlugin = PyRATBridge()
        self.mainPlugin.initGui()


def _writeLDPRELOAD():
    try:
        with open(os.path.join(
                os.environ["CONDA_PREFIX"], "etc", "conda",
                "activate.d", "qgis-activate.sh"
                ), "a") as script:
            script.write("export LD_PRELOAD=\"" + os.path.join(
                    os.environ["CONDA_PREFIX"], "lib",
                    "libedit.so") + " $LD_PRELOAD\"")
        with open(os.path.join(
                os.environ["CONDA_PREFIX"], "etc", "conda",
                "deactivate.d", "qgis-deactivate.sh"
                ), "a") as script:
            script.write("unset LD_PRELOAD")
        QMessageBox.warning(
                iface.mainWindow(),
                "Finish installation",
                "To finish the PyRATBridge installation please "
                "restart QGIS and reload your anaconda "
                "environment.")
    except Exception:
        QMessageBox.warning(
                iface.mainWindow(),
                "Exception during LD_PRELOAD configuration",
                "The installation is uncomplete and might crash on loading. "
                "Please reload the plugin or restart QGIS.")
