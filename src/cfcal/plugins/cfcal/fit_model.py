from nexpy.gui.utils import get_mainwindow, report_error
from nexpy.gui.widgets import NXCheckBox, NXComboBox
from nexusformat.nexus import NeXusError

from .base import B_PARAMETERS, SYMMETRY_PARAMETERS, CFParametersDialog

KERNELS = ['Gaussian', 'ChopperResolution']
CF_PARAMETERS = B_PARAMETERS + ['Hz', 'Hx', 'T']


def show_dialog(parent=None):
    try:
        dialog = FitModelDialog()
        dialog.show()
    except NeXusError as error:
        report_error("Fitting CF Model", error)


def convolution_available():
    """Return True if the Fit Panel supports convolve(model, kernel)."""
    try:
        from nexpy.gui import fitdialogs
        return hasattr(fitdialogs, 'convolve')
    except ImportError:
        return False


class FitModelDialog(CFParametersDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        action_buttons = self.action_buttons(('Load', self.load_parameters),
                                             ('Fit', self.fit_model),
                                             ('Save', self.write_parameters))
        self.convolve_checkbox = NXCheckBox('Convolve with',
                                            slot=self.set_convolution)
        self.kernel_combo = NXComboBox(items=KERNELS)
        if not convolution_available():
            self.convolve_checkbox.setEnabled(False)
            self.convolve_checkbox.setToolTip(
                "The installed version of NeXpy does not support convolution")
        convolve_layout = self.make_layout(self.convolve_checkbox,
                                           self.kernel_combo, align='center')
        self.set_layout(self.parameters.grid(header=False),
                        action_buttons, convolve_layout,
                        self.close_buttons(close=True))
        self.set_title('Fitting CF Model')

        self.fit_tab = None
        self.set_symmetry()
        self.set_convolution()

    def set_convolution(self):
        self.kernel_combo.setEnabled(self.convolve_checkbox.isChecked())

    def fit_model(self):
        """Initialize a CF model in a new tab of the NeXpy Fit Panel."""
        try:
            self.get_cf()  # validates the table values
            mainwindow = get_mainwindow()
            panel = mainwindow.panels.get('Fit')
            labels = list(panel.tabs) if panel else []
            mainwindow.fit_data()
            panel = mainwindow.panels.get('Fit')
            new = [label for label in (panel.tabs if panel else [])
                   if label not in labels]
            if not new:
                raise NeXusError("Select an NXdata group that is not already "
                                 "being fitted")
            tab = panel.tabs[new[0]]
            try:
                self.initialize_tab(tab)
            except Exception:
                panel.remove(new[0])
                raise
            self.fit_tab = tab
            self.capture_apply(tab)
        except Exception as error:
            report_error("Fitting CF Model", error)

    def initialize_tab(self, tab):
        """Add the CF model, and optionally a convolution, to a FitTab."""
        convolve = self.convolve_checkbox.isChecked()
        kernel = self.kernel_combo.selected
        if 'CF' not in tab.all_models:
            raise NeXusError("The CF model is not available in the Fit Panel")
        if convolve and kernel not in tab.all_models:
            raise NeXusError(f"The {kernel} model is not available")

        tab.model_combo.select('CF')
        tab.choose_model()
        tab.form_combo.select(self.parameters['RE'].value)
        tab.add_model()
        cf_model = tab.models[-1]
        prefix = cf_model['model'].prefix
        allowed = SYMMETRY_PARAMETERS[self.parameters['symmetry'].value]
        for name in CF_PARAMETERS:
            p = cf_model['parameters'][prefix + name]
            p.value = self.parameters[name].value
            if name in B_PARAMETERS:
                if name in allowed:
                    p.vary = True
                else:
                    p.value = 0.0
                    p.vary = False

        if convolve:
            # A line spectrum is convolved analytically with the kernel
            sigma = cf_model['parameters'][prefix + 'sigma']
            sigma.value = 0.0
            sigma.vary = False
            tab.model_combo.select(kernel)
            tab.choose_model()
            tab.add_model()
            expression = f"convolve({cf_model['name']}, {tab.models[-1]['name']})"
            tab.model = tab.eval_model(expression)
            tab.composite_model = expression
        tab.write_parameters()

    def capture_apply(self, tab):
        """Download the CF parameters when Apply is clicked in the tab."""
        apply = tab.apply

        def apply_tab():
            apply()
            try:
                self.fetch_parameters(tab)
            except Exception:
                pass

        tab.apply = apply_tab

    def fetch_parameters(self, tab):
        """Copy the CF parameters from a FitTab into the dialog table."""
        tab.read_parameters()
        prefix = tab.models[0]['model'].prefix
        parameters = tab.models[0]['parameters']
        for name in CF_PARAMETERS:
            value = parameters[prefix + name].value
            if value is not None and (name not in B_PARAMETERS
                                      or self.parameters[name].box.isEnabled()):
                self.parameters[name].value = value
