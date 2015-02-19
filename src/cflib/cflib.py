# CFlib  - Crystal Field Module
#
# Copyright (C) 2008-2015 R. Osborn, E. A. Goremychkin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
"""
CFlib: Python calculator for crystal fields

Crystal Field calculations using the Stevens Operator formalism.
"""

from scipy.linalg import eigh
import numpy as np

REs = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"]
Jf = [2.5, 4.0, 4.5, 4.0, 2.5, 0.0, 3.5, 6.0, 7.5, 8.0, 7.5, 6.0, 3.5]
gJ = [6.0/7.0, 4.0/5.0, 8.0/ 11.0, 3.0 / 5.0, 2.0 / 7.0, 0.0, 2.0, 3.0 / 2.0, 4.0 / 3.0,
      5.0 / 4.0, 6.0 / 5.0, 7.0 / 6.0, 8.0 / 7.0]
#<r^2> radial integral for 3+ f-electrons
r2 = [0.3666, 0.3380, 0.3120, 0.2917, 0.2728, 0.2569, 0.2428, 0.2302, 0.2188, 0.2085,
      0.1991, 0.1905, 0.1826]
#<r^4> radial integral for 3+ f-electrons
r4 = [0.3108, 0.2670, 0.2015, 0.1488, 0.1772, 0.1584, 0.1427, 0.1295, 0.1180, 0.1081,
      0.0996, 0.0921, 0.0854]
#<r^6> radial integral for 3+ f-electrons
r6 = [0.5119, 0.4150, 0.3300, 0.2787, 0.2317, 0.1985, 0.1720, 0.1505, 0.1328, 0.1810,
      0.1058, 0.0953, 0.0863]
#Second-degree Stevens factors for 3+ f-electrons
alphaJ = [-5.7143e-02, -2.1010e-02, -6.4279e-03, 7.7135e-03, 4.1270e-02, 0.0, 0.0,
          -1.0101e-02, -6.3492e-03, -2.2222e-03, 2.5397e-03, 1.0101e-02, 3.1746e-02]
#Fourth-degree Stevens factors for 3+ f-electrons
betaJ = [6.3492e-03, -7.3462e-04, -2.9111e-04, 4.0755e-04, 2.5012e-03, 0.0, 0.0,
         1.2244e-04, -5.9200e-05, -3.3300e-05, 4.4400e-05, 1.6325e-04, -1.7316e-03]
#Sixth-degree Stevens factors for 3+ f-electrons
gammaJ = [0.0, 6.0994e-05, -3.7988e-05, 6.6859e-04, 0.0, 0.0, 0.0, -1.1212e-06,
          1.0350e-06, -1.2937e-06, 2.0699e-06, -5.6061e-06, 1.4800e-04]


class CF(object):
    """
       Class defining the trivalent rare earth compound and its crystal field
       parameters
       
       There is only one basic object in CFlib, defining the trivalent rare earth,
       its crystal field parameters, and, if already diagonalized, the eigenvalues
       and eigenfunctions of the CF Hamiltonian.  
    """

    def __init__(self, RE=None, title=None, parfile=None):
        """Initializes the CF object or read it from a shelved file."""

        if RE:
            self.RE = RE
            try:
                self.index = REs.index(self.RE)
            except ValueError:
                print "CFerror: Rare Earth not recognized"
                return
        else:
            self.index = 1
        self.title = title
        self.Nf = self.index + 1
        self.J = Jf[self.index]
        self.size = int(2 * self.J + 1)
        self.gJ = gJ[self.index]
        self.r2 = r2[self.index]
        self.r4 = r4[self.index]
        self.r6 = r6[self.index]
        self.alphaJ = alphaJ[self.index]
        self.betaJ = betaJ[self.index]
        self.gammaJ = gammaJ[self.index]
        self.B20 = 0.0
        self.B40 = 0.0
        self.B60 = 0.0
        self.B22 = 0.0
        self.B42 = 0.0
        self.B62 = 0.0
        self.B43 = 0.0
        self.B63 = 0.0
        self.B44 = 0.0
        self.B64 = 0.0
        self.B66 = 0.0
        self.Hz = 0.0
        self.Hx = 0.0
        self.H = np.matrix(np.zeros((self.size, self.size), dtype='float64'))
        self.Jz = np.zeros((self.size, self.size), dtype='float64')
        self.Jp = np.zeros((self.size, self.size), dtype='float64')
        self.Jm = np.zeros((self.size, self.size), dtype='float64')
        self.TP = np.zeros((self.size, self.size), dtype='float64')
        self.EF = np.zeros((self.size, self.size), dtype='float64')
        self.EV = np.zeros(self.size, dtype='float64')
        self.peaks = []
        self.T = 0.0
        self.Tused = self.T
        self.resolution = 0.01
        self.threshold = 0.0001
        self.Jz_exp = 0.0
        self.Jx_exp = 0.0
        self.muz = 0.0
        self.mux = 0.0

        if parfile:
            import shelve

            try:
                db = shelve.open(parfile)
                if "CF" in db:
                    for item in db["CF"].__dict__.keys():
                        self.__dict__[item] = db["CF"].__dict__[item]
                db.close()
            except:
                print "CFerror: Error in reading file"

    def __str__(self):
        """Print the rare earth, CF parameters, and, if diagonalized, 
           the eigenvalues and eigenvectors."""

        if self.title:
            output = [self.title]
        else:
            output = []
        output.append("%s: Nf = %s, J = %s" % (self.RE, self.Nf, self.J))
        line = []
        if self.B20 != 0.0: line.append("B20 = %g" % self.B20)
        if self.B22 != 0.0: line.append("B22 = %g" % self.B22)
        if self.B40 != 0.0: line.append("B40 = %g" % self.B40)
        if self.B42 != 0.0: line.append("B42 = %g" % self.B42)
        if self.B43 != 0.0: line.append("B43 = %g" % self.B43)
        if self.B44 != 0.0: line.append("B44 = %g" % self.B44)
        if self.B60 != 0.0: line.append("B60 = %g" % self.B60)
        if self.B62 != 0.0: line.append("B62 = %g" % self.B62)
        if self.B63 != 0.0: line.append("B63 = %g" % self.B63)
        if self.B64 != 0.0: line.append("B64 = %g" % self.B64)
        if self.B66 != 0.0: line.append("B66 = %g" % self.B66)
        if self.Hz != 0.0: line.append("Hz = %g" % self.Hz)
        if self.Hx != 0.0: line.append("Hx = %g" % self.Hx)

        if line:
            output.append(" ".join(line))

        self.Ham()
        Jx_exp, Jz_exp, mux, muz = self.get_moments()

        if muz != 0.0 or mux != 0.0:
            output.append("<Jz> = %7.3f <muz> = %7.3f muB" % (Jz_exp, muz))
            output.append("<Jx> = %7.3f <mux> = %7.3f muB" % (Jx_exp, mux))

        if self.EV.any():
            output.append("Crystal Field Eigenvalues and Eigenfunctions")
            for i in range(self.EV.size):
                line = ["%8.3f:" % self.EV[i]]
                for j in range(self.EV.size):
                    if abs(self.EF[j, i]) > 0.0001:
                        Jz = j - self.J
                        if self.EF[j, i] < 0.0:
                            operator = "-"
                        else:
                            operator = "+"
                        line.append("%s%7.4f|%g>" % (operator, abs(self.EF[j, i]), Jz))
                output.append(" ".join(line))
        if self.peaks:
            output.append("Crystal Field Transitions")
            if self.T != self.Tused:
                self.get_peaks()
            output.append("Temperature: %g K" % self.T)
            for peak in self.peaks:
                output.append("Energy: %8.3f meV  Intensity: %8.4f" % peak)
        return "\n".join(output)

    def save(self):
        """Store the current object for later use."""

        import shelve

        db = None
        try:
            if self.title:
                db = shelve.open("%s.db" % self.title)
            else:
                db = shelve.open("CF.db")
            db["CF"] = self
        finally:
            if db:
                db.close()

    def read(self, parfile):
        """Load a saved version of the current object."""

        import shelve

        try:
            db = shelve.open(parfile)
            if "CF" in db:
                for item in db["CF"].__dict__.keys():
                    self.__dict__[item] = db["CF"].__dict__[item]
            db.close()
        except:
            print "CFerror: Error in reading file"

    def CFham(self):
        """Determine the CF Hamiltonian based on the input parameters."""

        J = self.J
        H = np.matrix(np.zeros((self.size, self.size), dtype='float64'))

        for m in range(self.size):
            mJ = m - J
            O20 = 3 * mJ ** 2 - J * (J + 1)
            O40 = 35 * mJ ** 4 - 30 * J * (J + 1) * mJ ** 2 + 25 * mJ ** 2 - 6 * J * (J + 1) + 3 * (J * (J + 1)) ** 2
            O60 = 231 * mJ ** 6 - 315 * J * (J + 1) * mJ ** 4 + 735 * mJ ** 4 + 105 * (J * (J + 1) * mJ) ** 2 \
                  - 525 * J * (J + 1) * mJ ** 2 + 294 * mJ ** 2 - 5 * (J * (J + 1)) ** 3 + 40 * (J * (J + 1)) ** 2 - 60 * J * (J + 1)
            H[m, m] += self.B20 * O20 + self.B40 * O40 + self.B60 * O60

        for m in range(self.size - 2):
            mJ = m - J
            n = m + 2
            nJ = mJ + 2
            O22 = 0.5 * np.sqrt((J * (J + 1) - nJ * (nJ - 1)) * (J * (J + 1) - (nJ - 1) * (nJ - 2)))
            O42 = (3.5 * (mJ ** 2 + nJ ** 2) - J * (J + 1) - 5) * O22
            O62 = (16.5 * (mJ ** 4 + nJ ** 4) - 9 * (mJ ** 2 + nJ ** 2) * J * (J + 1) - 61.5 * (mJ ** 2 + nJ ** 2)
                   + (J * (J + 1)) ** 2 + 10 * J * (J + 1) + 102) * O22
            H[m, n] += self.B22 * O22 + self.B42 * O42 + self.B62 * O62
            H[n, m] = H[m, n]

        for m in range(self.size - 3):
            mJ = m - J
            n = m + 3
            nJ = mJ + 3
            A = 1.0
            for k in range(3):
                A *= J * (J + 1) - (nJ - k) * (nJ - k - 1)
            O43 = 0.25 * np.sqrt(A) * (mJ + nJ)
            O63 = 0.25 * (11 * (mJ ** 3 + nJ ** 3) - 3 * (mJ + nJ) * J * (J + 1) - 59 * (mJ + nJ)) * np.sqrt(A)
            H[m, n] += self.B43 * O43 + self.B63 * O63
            H[n, m] = H[m, n]

        for m in range(self.size - 4):
            mJ = m - J
            n = m + 4
            nJ = mJ + 4
            O44 = 1.0
            for k in range(4):
                O44 *= J * (J + 1) - (nJ - k) * (nJ - k - 1)
            O44 = 0.5 * np.sqrt(O44)
            O64 = (5.5 * (mJ ** 2 + nJ ** 2) - J * (J + 1) - 38) * O44
            H[m, n] += self.B44 * O44 + self.B64 * O64
            H[n, m] = H[m, n]

        for m in range(self.size - 6):
            mJ = m - J
            n = m + 6
            nJ = mJ + 6
            O66 = 1.0
            for k in range(6):
                O66 *= J * (J + 1) - (nJ - k) * (nJ - k - 1)
            O66 = 0.5 * np.sqrt(O66)
            H[m, n] += self.B66 * O66
            H[n, m] = H[m, n]

        return H

    def MFham(self):
        """Determine the Zeeman terms to the Hamiltonian."""

        J = self.J
        H = np.matrix(np.zeros((self.size, self.size), dtype='float64'))

        for m in range(self.size):
            mJ = m - J
            H[m, m] -= self.Hz * mJ

        for m in range(self.size - 1):
            mJ = m - J
            n = m + 1
            nJ = mJ + 1
            H[m, n] -= 0.5 * self.Hx * np.sqrt((J * (J + 1) - mJ * nJ))
            H[n, m] = H[m, n]

        return H

    def Ham(self):
        """Returns the total Hamiltonian including CF and Zeeman terms."""

        return self.CFham() + self.MFham()

    def diagonalize(self):
        """Determine the eigenvalues and eigenfunctions of the total Hamiltonian."""

        H = self.Ham()
        self.EV, self.EF = eigh(H)
        self.EV = self.EV - self.EV.min()

    def TPS(self):
        """Determine the dipole matrix elements for the total Hamiltonian."""

        J = self.J

        self.Jz = np.zeros((self.size, self.size), dtype='float64')
        self.Jp = np.zeros((self.size, self.size), dtype='float64')
        self.Jm = np.zeros((self.size, self.size), dtype='float64')
        for m in range(self.size):
            for n in range(self.size):
                for k in range(self.size):
                    kJ = k - J
                    self.Jz[m, n] += self.EF[k, m] * self.EF[k, n] * kJ
                for k in range(self.size - 1):
                    kJ = k - J
                    l = k + 1
                    lJ = l - J
                    self.Jp[m, n] = self.Jp[m, n] + \
                                    self.EF[l, m] * self.EF[k, n] * np.sqrt(J * (J + 1) - kJ * lJ)
                for k in range(1, self.size):
                    kJ = k - J
                    l = k - 1
                    lJ = l - J
                    self.Jm[m, n] = self.Jm[m, n] + \
                                    self.EF[l, m] * self.EF[k, n] * np.sqrt(J * (J + 1) - kJ * lJ)
                self.TP[m, n] = (2 * self.Jz[m, n] ** 2 + self.Jp[m, n] ** 2 + self.Jm[m, n] ** 2) / 3

    def get_peaks(self, T=None, Hx=None, Hz=None):
        """Determine the peak intensities from the total Hamiltonian."""

        if Hx:
            old_Hx = self.Hx
            self.Hx = Hx
        if Hz:
            old_Hz = self.Hz
            self.Hz = Hz
        self.EFS()
        self.TPS()

        BF = np.zeros(self.size, dtype='float64')

        if T is None: T = self.T
        kT = T / 11.6045
        if kT <= 0.0:
            BF[0] = 1.0
        else:
            Z = sum(np.exp(-self.EV / kT))
            BF = np.exp(-self.EV / kT) / Z

        peaks = []
        for n in range(self.size):
            for m in range(self.size):
                ET = self.EV[m] - self.EV[n]
                IT = self.TP[m, n] * BF[n]
                if IT > 0.0:
                    peaks.append((ET, IT))

        self.peaks = []
        for k in range(len(peaks)):
            ETk, ITk = peaks[k]
            if ITk > 0.0:
                sum_peaks = ETk * ITk
                for l in range(len(peaks)):
                    if k != l:
                        ETl, ITl = peaks[l]
                        if ITl > 0.0 and abs(ETk - ETl) < self.resolution:
                            ITk = ITk + ITl
                            sum_peaks = sum_peaks + ETl * ITl
                            peaks[l] = (ETl, 0.0)
                if ITk > self.threshold:
                    ETk = sum_peaks / ITk
                    self.peaks.append((ETk, ITk))
        self.peaks.sort()
        self.Tused = T
        if Hx:
            self.Hx = old_Hx
        if Hz:
            self.Hz = old_Hz

        return self.peaks

    def spectrum(self, eps=None, sigma=None, T=None, Hx=None, Hz=None):
        """Calculates the neutron scattering cross section.
        """

        if T is None:
            T = self.T

        peaks = self.get_peaks(T, Hx, Hz)

        if eps is None: eps = np.linspace(-1.1 * self.EV[-1], 1.1 * self.EV[-1], 501)
        if sigma is None: sigma = 0.01 * (max(eps) - min(eps))

        S = np.zeros(eps.size, dtype='float64')

        for EV, IT in peaks:
            S += IT * np.exp(-(eps - EV) ** 2 / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))
        S *= 72.65 * self.gJ ** 2

        return S

    def NXspectrum(self, eps=None, sigma=None, T=None, Hx=None, Hz=None):
        """Returns the neutron scattering cross section as a NXentry group"""

        from nexusformat.nexus import NXfield, NXentry, NXsample, NXdata

        S = self.spectrum(eps, sigma, T, Hx, Hz)
        entry = NXentry()
        entry.title = "Crystal Field Spectra for %s at %s K" % (self.title, T)
        entry.sample = NXsample()
        entry.sample.temperature = T
        entry.sample.temperature.units = "K"
        entry.data = NXdata(NXfield(S, name="intensity", units="mb/sr/meV"),
                            NXfield(eps, name="energy_transfer", units="meV"))
        return entry

    def get_moments(self, T=None):
        """Calculate the magnetic moments of the CF model."""

        if T is None:
            T = self.T
        if T != self.Tused:
            self.get_peaks(T)
        kT = T / 11.6045
        if kT > 0.0:
            Jz_exp = 0.0
            Jx_exp = 0.0
            Z = sum(np.exp(-self.EV / kT))
            for i in range(self.EV.size):
                Jz_exp = Jz_exp + self.Jz[i, i] * np.exp(-self.EV[i] / kT)
                Jx_exp += 0.5 * (self.Jp[i, i] + self.Jm[i, i]) * np.exp(-self.EV[i] / kT)
            Jz_exp = Jz_exp / Z
            Jx_exp = Jx_exp / Z
        else:
            Jz_exp = sum(self.Jz[self.EV == 0.0, self.EV == 0.0]) / self.EV[self.EV == 0.0].size
            Jx_exp = sum(0.5 * (self.Jp[self.EV == 0.0, self.EV == 0.0]
                                + self.Jm[self.EV == 0.0, self.EV == 0.0])) \
                     / self.EV[self.EV == 0.0].size
        self.Jz_exp = Jz_exp
        self.muz = self.gJ * self.Jz_exp
        self.Jx_exp = Jx_exp
        self.mux = self.gJ * self.Jx_exp

        return self.Jx_exp, self.Jz_exp, self.mux, self.muz

    def chi(self, T=None):
        """Calculate the susceptibility at the specified temperatures."""

        if T is None:
            T = np.linspace(1.0, 300.0, 300)
        kT = T / 11.6045
        Z = sum(np.exp(-self.EV / kT))
        ChiC_zz = 0.0 * T
        ChiC_xx = 0.0 * T
        ChiV_zz = 0.0 * T
        ChiV_xx = 0.0 * T
        for m in range(self.EV.size):
            for n in range(self.EV.size):
                if all(abs(self.EV[n] - self.EV[m]) < 0.00001 * kT):
                    ChiC_zz += self.Jz[m, n] ** 2 * np.exp(-self.EV[m] / kT)
                    ChiC_xx += 0.25 * (self.Jp[m, n] ** 2 + self.Jm[m, n] ** 2) \
                               * np.exp(-self.EV[m] / kT)
                else:
                    ChiV_zz += 2 * self.Jz[m, n] ** 2 * np.exp(-self.EV[m] / kT) \
                               / (self.EV[n] - self.EV[m])
                    ChiV_xx += 0.5 * (self.Jp[m, n] ** 2 + self.Jm[m, n] ** 2) \
                               * np.exp(-self.EV[m] / kT) / (self.EV[n] - self.EV[m])
        ChiC_zz = self.gJ ** 2 * ChiC_zz / (kT * Z)
        ChiC_xx = self.gJ ** 2 * ChiC_xx / (kT * Z)
        ChiV_zz = self.gJ ** 2 * ChiV_zz / Z
        ChiV_xx = self.gJ ** 2 * ChiV_xx / Z

        return ChiC_zz, ChiC_xx, ChiV_zz, ChiV_xx

    def NXchi(self, T=None):
        """Returns the susceptibility as a NeXus NXentry"""

        from nexusformat.nexus import NXfield, NXentry, NXdata

        ChiC_zz, ChiC_xx, ChiV_zz, ChiV_xx = self.chi(T)

        entry = NXentry()
        entry.title = self.title
        temperature = NXfield(T, name="temperature")
        temperature.units = "K"

        chi = NXfield(ChiC_zz + ChiV_zz + 2 * (ChiC_xx + ChiV_xx), name="chi")

        entry.chi = NXdata(chi, temperature)
        entry.chi.title = "Susceptibility of %s" % self.title
        entry.invchi = NXdata(1 / chi, temperature)
        entry.invchi.title = "Inverse Susceptibility of %s" % self.title
        entry.chiz = NXdata(ChiC_zz + ChiV_zz, temperature)
        entry.chiz.title = "Susceptibility of %s (z-axis)" % self.title
        entry.chix = NXdata(ChiC_xx + ChiV_xx, temperature)
        entry.chix.title = "Susceptibility of %s (x-axis)" % self.title

        return entry

