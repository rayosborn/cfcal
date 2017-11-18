import numpy as np
from nexpy.gui.datadialogs import BaseDialog, GridParameters
from nexpy.gui.plotview import plotview
from nexpy.gui.utils import report_error
from nexusformat.nexus import *
from cfcal.cfcal import CF

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

