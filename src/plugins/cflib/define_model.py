import numpy as np
from nexpy.gui.datadialogs import BaseDialog, GridParameters
from nexpy.gui.plotview import plotview
from nexpy.gui.utils import report_error
from nexusformat.nexus import *
from cflib.cflib import CF

def show_dialog(parent=None):
    try:
        dialog = DefineModelDialog()
        dialog.show()
    except NeXusError as error:
        report_error("Defining CF Model", error)


class DefineModelDialog(BaseDialog):

    def __init__(self, parent=None):
        super(DefineModelDialog, self).__init__(parent)

        node = self.get_node()
        self.root = node.nxroot
        
        symmetries = ['cubic', 'tetragonal', 'orthorhombic', 'hexagonal', 
                      'monoclinic', 'triclinic']

        rare_earths = ['Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 
                       'Ho', 'Er', 'Tm', 'Yb']

        self.parameters = GridParameters()
        self.parameters.add('symmetry', symmetries, 'Symmetry') 

        action_buttons = self.action_buttons(('Plot', self.plot_lattice),
                                             ('Save', self.write_parameters))
        self.set_layout(self.entry_layout, self.parameters.grid(), 
                        action_buttons, self.close_buttons())
        self.set_title('Defining CF Model')


    def cf_grid(self):
        parameters = []
        if self.symmetry == 'cubic':
            parameters
        
        self.B20_box = QtGui.QLineEdit()
        self.B22_box = QtGui.QLineEdit()
        self.B40_box = QtGui.QLineEdit()
        self.B42_box = QtGui.QLineEdit()
        self.B43_box = QtGui.QLineEdit()
        self.B44_box = QtGui.QLineEdit()
        self.B60_box = QtGui.QLineEdit()
        self.B62_box = QtGui.QLineEdit()
        self.B63_box = QtGui.QLineEdit()
        self.B64_box = QtGui.QLineEdit()
        self.B66_box = QtGui.QLineEdit()
        self.Hz_box = QtGui.QLineEdit()
        self.Hx_box = QtGui.QLineEdit()
        grid.addWidget(QtGui.QLabel('B20:'), 0, 0)
        grid.addWidget(QtGui.QLabel('B22:'), 0, 0)
        grid.addWidget(QtGui.QLabel('Unit Cell - a (Ang):'), 1, 0)
        grid.addWidget(QtGui.QLabel('Unit Cell - b (Ang):'), 2, 0)
        grid.addWidget(QtGui.QLabel('Unit Cell - c (Ang):'), 3, 0)
        grid.addWidget(QtGui.QLabel('Unit Cell - alpha (deg):'), 4, 0)
        grid.addWidget(QtGui.QLabel('Unit Cell - beta (deg):'), 5, 0)
        grid.addWidget(QtGui.QLabel('Unit Cell - gamma (deg):'), 6, 0)
        grid.addWidget(QtGui.QLabel('Wavelength (Ang):'), 7, 0)
        grid.addWidget(QtGui.QLabel('Distance (mm):'), 8, 0)
        grid.addWidget(QtGui.QLabel('Yaw (deg):'), 9, 0)
        grid.addWidget(QtGui.QLabel('Pitch (deg):'), 10, 0)
        grid.addWidget(QtGui.QLabel('Roll (deg):'), 11, 0)
 

    def update_parameter(self, box, value):
        if value is not None:
            box.setText(str(value))

    def update_parameters(self):
        self.update_parameter(self.unitcell_a_box, self.refine.a)
        self.update_parameter(self.unitcell_b_box, self.refine.b)
        self.update_parameter(self.unitcell_c_box, self.refine.c)
        self.update_parameter(self.unitcell_alpha_box, self.refine.alpha)
        self.update_parameter(self.unitcell_beta_box, self.refine.beta)
        self.update_parameter(self.unitcell_gamma_box, self.refine.gamma)
        self.update_parameter(self.wavelength_box, self.refine.wavelength)
        self.update_parameter(self.distance_box, self.refine.distance)
        self.update_parameter(self.yaw_box, self.refine.yaw)
        self.update_parameter(self.pitch_box, self.refine.pitch)
        self.update_parameter(self.roll_box, self.refine.roll)
        self.update_parameter(self.xc_box, self.refine.xc)
        self.update_parameter(self.yc_box, self.refine.yc)

    @property
    def symmetry(self):
        return self.symmetry_box.currentText()

    def set_symmetry(self):
        self.refine.symmetry = self.get_symmetry()
        self.refine.set_symmetry()
        self.update_parameters()
        if self.refine.symmetry == 'cubic':
            self.unitcell_b_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_c_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_alpha_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_beta_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_gamma_checkbox.setCheckState(QtCore.Qt.Unchecked)
        elif self.refine.symmetry == 'tetragonal':
            self.unitcell_b_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_alpha_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_beta_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_gamma_checkbox.setCheckState(QtCore.Qt.Unchecked)
        elif self.refine.symmetry == 'orthorhombic':
            self.unitcell_alpha_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_beta_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_gamma_checkbox.setCheckState(QtCore.Qt.Unchecked)
        elif self.refine.symmetry == 'hexagonal':
            self.unitcell_b_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_alpha_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_beta_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_gamma_checkbox.setCheckState(QtCore.Qt.Unchecked)
        elif self.refine.symmetry == 'monoclinic':
            self.unitcell_alpha_checkbox.setCheckState(QtCore.Qt.Unchecked)
            self.unitcell_gamma_checkbox.setCheckState(QtCore.Qt.Unchecked)

    def guess_symmetry(self):
        self.refine.a, self.refine.b, self.refine.c, \
            self.refine.alpha, self.refine.beta, self.refine.gamma = self.get_lattice_parameters()
        return self.refine.guess_symmetry()

    def get_lattice_parameters(self):
        return (np.float32(self.unitcell_a_box.text()),
                np.float32(self.unitcell_b_box.text()),
                np.float32(self.unitcell_c_box.text()),
                np.float32(self.unitcell_alpha_box.text()),
                np.float32(self.unitcell_beta_box.text()),
                np.float32(self.unitcell_gamma_box.text()))

    def get_wavelength(self):
        return np.float32(self.wavelength_box.text())

    def get_distance(self):
        return np.float32(self.distance_box.text())

    def get_tilts(self):
        return (np.float32(self.yaw_box.text()),
                np.float32(self.pitch_box.text()),
                np.float32(self.roll_box.text()))

    def get_centers(self):
        return np.float32(self.xc_box.text()), np.float32(self.yc_box.text())

    def get_polar_max(self):
        return np.float32(self.polar_box.text())

    def set_polar_max(self):
        self.refine.set_polar_max(self.get_polar_max())

    def initialize_fit(self):
        self.refine.parameters = []
        if self.unitcell_a_checkbox.isChecked():
            self.refine.parameters.append({'a':self.refine.a})
        if self.unitcell_b_checkbox.isChecked():
            self.refine.parameters.append({'b':self.refine.b})
        if self.unitcell_c_checkbox.isChecked():
            self.refine.parameters.append({'c':self.refine.c})
        if self.unitcell_alpha_checkbox.isChecked():
            self.refine.parameters.append({'alpha':self.refine.alpha})
        if self.unitcell_beta_checkbox.isChecked():
            self.refine.parameters.append({'beta':self.refine.beta})
        if self.unitcell_gamma_checkbox.isChecked():
            self.refine.parameters.append({'gamma':self.refine.gamma})
        if self.wavelength_checkbox.isChecked():
            self.refine.parameters.append({'wavelength':self.refine.wavelength})
        if self.distance_checkbox.isChecked():
            self.refine.parameters.append({'distance':self.refine.distance})
        if self.yaw_checkbox.isChecked():
            self.refine.parameters.append({'yaw':self.refine.yaw})
        if self.pitch_checkbox.isChecked():
            self.refine.parameters.append({'pitch':self.refine.pitch})
        if self.roll_checkbox.isChecked():
            self.refine.parameters.append({'roll':self.refine.roll})
        if self.xc_checkbox.isChecked():
            self.refine.parameters.append({'xc':self.refine.xc})
        if self.yc_checkbox.isChecked():
            self.refine.parameters.append({'yc':self.refine.yc})

    def plot_peaks(self):
        self.refine.polar_max = self.get_polar_max()
        self.refine.plot_peaks(self.refine.x, self.refine.y)
        self.refine.plot_rings()

    def refine_parameters(self):
        self.initialize_fit()
        self.refine.refine_parameters()
        self.update_parameters()

    def write_parameters(self):
        try:
            polar_angles, azimuthal_angles = self.refine.calculate_angles(
                                                 self.refine.xp, self.refine.yp)
            self.refine.write_angles(polar_angles, azimuthal_angles)
            self.refine.write_parameters()
        except NeXusError as error:
            report_error('Refining Lattice', error)
