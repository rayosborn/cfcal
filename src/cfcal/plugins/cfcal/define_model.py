from nexpy.gui.plotview import NXPlotView, plotviews
from nexpy.gui.pyqt import QtCore
from nexpy.gui.utils import report_error
from nexpy.gui.widgets import NXCheckBox
from nexusformat.nexus import NeXusError

from .base import CFParametersDialog

# Each quantity maps to the plot window label, the CF method returning an
# NXentry, and the name of the NXdata group within that entry to plot.
QUANTITIES = {
    'Spectrum': ('Crystal Field Spectrum', 'NXspectrum', 'data'),
    'Susceptibility': ('Crystal Field Susceptibility', 'NXchi', 'chi'),
    'Inverse Susceptibility': ('Crystal Field Inverse Susceptibility',
                               'NXchi', 'invchi'),
    'Specific Heat': ('Crystal Field Specific Heat', 'NXthermo',
                      'specific_heat'),
    'Entropy': ('Crystal Field Entropy', 'NXthermo', 'entropy'),
}


def show_dialog(parent=None):
    try:
        dialog = DefineModelDialog()
        dialog.show()
    except NeXusError as error:
        report_error("Defining CF Model", error)


class DefineModelDialog(CFParametersDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.parameters.add('quantity', list(QUANTITIES), 'Plot',
                            slot=self.update_overplot_checkbox)

        action_buttons = self.action_buttons(('Load', self.load_parameters),
                                             ('Plot', self.plot_quantity),
                                             ('Save', self.write_parameters))
        self.overplot_checkbox = NXCheckBox('Overplot')
        self.line_checkbox = NXCheckBox('Plot as Line')
        checkbox_layout = self.make_layout(self.overplot_checkbox,
                                           self.line_checkbox, align='center')
        self.set_layout(self.parameters.grid(header=False),
                        action_buttons, checkbox_layout,
                        self.close_buttons(close=True))
        self.set_title('Defining CF Model')

        self.set_symmetry()
        self.update_overplot_checkbox()

    @property
    def plot_label(self):
        return QUANTITIES[self.parameters['quantity'].value][0]

    def update_overplot_checkbox(self):
        exists = self.plot_label in plotviews
        self.overplot_checkbox.setEnabled(exists)
        if not exists:
            self.overplot_checkbox.setChecked(False)

    def changeEvent(self, event):
        super().changeEvent(event)
        if (event.type() == QtCore.QEvent.ActivationChange
                and self.isActiveWindow()):
            self.update_overplot_checkbox()

    def plot_quantity(self):
        quantity = self.parameters['quantity'].value
        label, method, group = QUANTITIES[quantity]
        try:
            cf = self.get_cf()
            entry = getattr(cf, method)()
            opts = {}
            if (self.overplot_checkbox.isEnabled()
                    and self.overplot_checkbox.isChecked()):
                opts['over'] = True
            if self.line_checkbox.isChecked():
                opts['marker'] = 'None'
                opts['linestyle'] = '-'
            if label in plotviews:
                plotview = plotviews[label]
            else:
                plotview = NXPlotView(label)
            plotview.plot(entry[group], **opts)
            self.update_overplot_checkbox()
        except Exception as error:
            report_error(f"Plotting CF {quantity}", error)
