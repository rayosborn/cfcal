from nexpy.gui.dialogs import GridParameters, NXDialog
from nexpy.gui.plotview import NXPlotView, plotviews
from nexpy.gui.pyqt import QtCore, getOpenFileName, getSaveFileName
from nexpy.gui.utils import report_error
from nexpy.gui.widgets import NXCheckBox
from nexusformat.nexus import NeXusError

PLOT_LABEL = 'Crystal Field Spectrum'

B_PARAMETERS = ['B20', 'B22', 'B40', 'B42', 'B43', 'B44',
                'B60', 'B62', 'B63', 'B64', 'B66']

SYMMETRY_PARAMETERS = {
    'cubic': ['B40', 'B44', 'B60', 'B64'],
    'tetragonal': ['B20', 'B40', 'B44', 'B60', 'B64'],
    'orthorhombic': ['B20', 'B22', 'B40', 'B42', 'B44',
                     'B60', 'B62', 'B64', 'B66'],
    'hexagonal': ['B20', 'B40', 'B60', 'B66'],
    'monoclinic': ['B20', 'B22', 'B40', 'B42', 'B44',
                   'B60', 'B62', 'B64', 'B66'],
    'triclinic': B_PARAMETERS,
}


def show_dialog(parent=None):
    try:
        dialog = DefineModelDialog()
        dialog.show()
    except NeXusError as error:
        report_error("Defining CF Model", error)


class DefineModelDialog(NXDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        symmetries = ['cubic', 'tetragonal', 'orthorhombic', 'hexagonal',
                      'monoclinic', 'triclinic']

        self.rare_earths = ['Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                            'Dy', 'Ho', 'Er', 'Tm', 'Yb']

        self.parameters = GridParameters()
        self.parameters.add('RE', self.rare_earths, 'Rare Earth')
        self.parameters.add('name', '', 'Name')
        self.parameters.add('symmetry', symmetries, 'Symmetry',
                            slot=self.set_symmetry)
        for name in B_PARAMETERS:
            self.parameters.add(name, 0.0, name)
        self.parameters.add('Hz', 0.0, 'Hz')
        self.parameters.add('Hx', 0.0, 'Hx')
        self.parameters.add('T', 0.0, 'Temperature (K)')

        action_buttons = self.action_buttons(('Load', self.load_parameters),
                                             ('Plot', self.plot_spectrum),
                                             ('Save', self.write_parameters))
        self.overplot_checkbox = NXCheckBox('Overplot Spectrum')
        self.line_checkbox = NXCheckBox('Plot as Line')
        checkbox_layout = self.make_layout(self.overplot_checkbox,
                                           self.line_checkbox, align='center')
        self.set_layout(self.parameters.grid(header=False),
                        action_buttons, checkbox_layout,
                        self.close_buttons(close=True))
        self.set_title('Defining CF Model')

        self.set_symmetry()
        self.update_overplot_checkbox()

    def set_symmetry(self):
        enabled = SYMMETRY_PARAMETERS[self.parameters['symmetry'].value]
        for name in B_PARAMETERS:
            self.parameters[name].box.setEnabled(name in enabled)

    def update_overplot_checkbox(self):
        exists = PLOT_LABEL in plotviews
        self.overplot_checkbox.setEnabled(exists)
        if not exists:
            self.overplot_checkbox.setChecked(False)

    def changeEvent(self, event):
        super().changeEvent(event)
        if (event.type() == QtCore.QEvent.ActivationChange
                and self.isActiveWindow()):
            self.update_overplot_checkbox()

    def infer_symmetry(self, nonzero):
        """Return the smallest symmetry whose parameters include `nonzero`."""
        for symmetry in sorted(SYMMETRY_PARAMETERS,
                               key=lambda s: len(SYMMETRY_PARAMETERS[s])):
            if set(nonzero).issubset(SYMMETRY_PARAMETERS[symmetry]):
                return symmetry
        return 'triclinic'

    def set_cf(self, cf):
        self.parameters['RE'].value = cf.RE
        self.parameters['name'].value = cf.name or ''
        nonzero = [name for name in B_PARAMETERS
                  if getattr(cf, name) != 0.0]
        self.parameters['symmetry'].value = self.infer_symmetry(nonzero)
        self.set_symmetry()
        for name in B_PARAMETERS:
            self.parameters[name].value = getattr(cf, name)
        self.parameters['Hz'].value = cf.Hz
        self.parameters['Hx'].value = cf.Hx
        self.parameters['T'].value = cf.T

    def get_cf(self):
        from cfcal import CF

        cf = CF(RE=self.parameters['RE'].value,
                name=self.parameters['name'].value or None)
        for name in B_PARAMETERS:
            if self.parameters[name].box.isEnabled():
                setattr(cf, name, self.parameters[name].value)
            else:
                setattr(cf, name, 0.0)
        cf.Hz = self.parameters['Hz'].value
        cf.Hx = self.parameters['Hx'].value
        cf.T = self.parameters['T'].value
        return cf

    def load_parameters(self):
        from cfcal import CF

        try:
            fname = getOpenFileName(self, "Choose a Filename")
            if fname:
                self.set_cf(CF(parfile=fname))
        except Exception as error:
            report_error("Loading CF Parameters", error)

    def plot_spectrum(self):
        try:
            cf = self.get_cf()
            entry = cf.NXspectrum()
            opts = {}
            if (self.overplot_checkbox.isEnabled()
                    and self.overplot_checkbox.isChecked()):
                opts['over'] = True
            if self.line_checkbox.isChecked():
                opts['marker'] = 'None'
                opts['linestyle'] = '-'
            if PLOT_LABEL in plotviews:
                plotview = plotviews[PLOT_LABEL]
            else:
                plotview = NXPlotView(PLOT_LABEL)
            plotview.plot(entry.data, **opts)
            self.update_overplot_checkbox()
        except Exception as error:
            report_error("Plotting CF Spectrum", error)

    def write_parameters(self):
        try:
            cf = self.get_cf()
            fname = getSaveFileName(self, "Choose a Filename",
                                    f"{cf.name or 'model'}.cfg")
            if fname:
                cf.save(fname)
        except Exception as error:
            report_error("Saving CF Parameters", error)
