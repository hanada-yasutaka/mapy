# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'energydomain_dialog.ui'
#
# Created: Fri Jul 25 11:53:32 2014
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_EnergyDomainDialog(object):
    def setupUi(self, EnergyDomainDialog):
        EnergyDomainDialog.setObjectName(_fromUtf8("EnergyDomainDialog"))
        EnergyDomainDialog.resize(440, 348)
        self.CloseButton = QtGui.QPushButton(EnergyDomainDialog)
        self.CloseButton.setGeometry(QtCore.QRect(320, 310, 114, 32))
        self.CloseButton.setAutoDefault(False)
        self.CloseButton.setObjectName(_fromUtf8("CloseButton"))
        self.groupBox_System = QtGui.QGroupBox(EnergyDomainDialog)
        self.groupBox_System.setGeometry(QtCore.QRect(10, 30, 191, 151))
        self.groupBox_System.setObjectName(_fromUtf8("groupBox_System"))
        self.GetEigenButton = QtGui.QPushButton(self.groupBox_System)
        self.GetEigenButton.setGeometry(QtCore.QRect(100, 120, 91, 32))
        self.GetEigenButton.setAutoDefault(False)
        self.GetEigenButton.setObjectName(_fromUtf8("GetEigenButton"))
        self.lineEdit_dim = QtGui.QLineEdit(self.groupBox_System)
        self.lineEdit_dim.setGeometry(QtCore.QRect(70, 32, 41, 21))
        self.lineEdit_dim.setObjectName(_fromUtf8("lineEdit_dim"))
        self.label = QtGui.QLabel(self.groupBox_System)
        self.label.setGeometry(QtCore.QRect(10, 32, 60, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.lineEdit_qdomain = QtGui.QLineEdit(self.groupBox_System)
        self.lineEdit_qdomain.setGeometry(QtCore.QRect(40, 64, 71, 21))
        self.lineEdit_qdomain.setObjectName(_fromUtf8("lineEdit_qdomain"))
        self.label_3 = QtGui.QLabel(self.groupBox_System)
        self.label_3.setGeometry(QtCore.QRect(10, 96, 16, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_2 = QtGui.QLabel(self.groupBox_System)
        self.label_2.setGeometry(QtCore.QRect(10, 64, 16, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.lineEdit_pdomain = QtGui.QLineEdit(self.groupBox_System)
        self.lineEdit_pdomain.setGeometry(QtCore.QRect(40, 96, 72, 21))
        self.lineEdit_pdomain.setObjectName(_fromUtf8("lineEdit_pdomain"))
        self.comboBox_SelectAlgorism = QtGui.QComboBox(EnergyDomainDialog)
        self.comboBox_SelectAlgorism.setGeometry(QtCore.QRect(10, 10, 381, 26))
        self.comboBox_SelectAlgorism.setObjectName(_fromUtf8("comboBox_SelectAlgorism"))
        self.comboBox_SelectAlgorism.addItem(_fromUtf8(""))
        self.comboBox_SelectAlgorism.addItem(_fromUtf8(""))
        self.comboBox_SelectAlgorism.addItem(_fromUtf8(""))
        self.groupBox_Hsm = QtGui.QGroupBox(EnergyDomainDialog)
        self.groupBox_Hsm.setEnabled(False)
        self.groupBox_Hsm.setGeometry(QtCore.QRect(210, 30, 221, 271))
        self.groupBox_Hsm.setObjectName(_fromUtf8("groupBox_Hsm"))
        self.lineEdit_hsmqrange = QtGui.QLineEdit(self.groupBox_Hsm)
        self.lineEdit_hsmqrange.setGeometry(QtCore.QRect(50, 60, 71, 21))
        self.lineEdit_hsmqrange.setObjectName(_fromUtf8("lineEdit_hsmqrange"))
        self.lineEdit_hsmprange = QtGui.QLineEdit(self.groupBox_Hsm)
        self.lineEdit_hsmprange.setGeometry(QtCore.QRect(50, 90, 71, 21))
        self.lineEdit_hsmprange.setObjectName(_fromUtf8("lineEdit_hsmprange"))
        self.label_4 = QtGui.QLabel(self.groupBox_Hsm)
        self.label_4.setGeometry(QtCore.QRect(20, 60, 20, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_6 = QtGui.QLabel(self.groupBox_Hsm)
        self.label_6.setGeometry(QtCore.QRect(20, 90, 21, 16))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.label_7 = QtGui.QLabel(self.groupBox_Hsm)
        self.label_7.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.lineEdit_hsmgrid = QtGui.QLineEdit(self.groupBox_Hsm)
        self.lineEdit_hsmgrid.setGeometry(QtCore.QRect(120, 30, 61, 21))
        self.lineEdit_hsmgrid.setObjectName(_fromUtf8("lineEdit_hsmgrid"))
        self.drawButton = QtGui.QPushButton(self.groupBox_Hsm)
        self.drawButton.setGeometry(QtCore.QRect(110, 120, 114, 32))
        self.drawButton.setAutoDefault(False)
        self.drawButton.setObjectName(_fromUtf8("drawButton"))
        self.label_5 = QtGui.QLabel(self.groupBox_Hsm)
        self.label_5.setGeometry(QtCore.QRect(29, 121, 21, 16))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.lineEdit_QuantumNumber = QtGui.QLineEdit(self.groupBox_Hsm)
        self.lineEdit_QuantumNumber.setGeometry(QtCore.QRect(50, 120, 61, 21))
        self.lineEdit_QuantumNumber.setObjectName(_fromUtf8("lineEdit_QuantumNumber"))
        self.label_9 = QtGui.QLabel(self.groupBox_Hsm)
        self.label_9.setGeometry(QtCore.QRect(10, 190, 101, 16))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.lineEdit_contourlinspace = QtGui.QLineEdit(self.groupBox_Hsm)
        self.lineEdit_contourlinspace.setGeometry(QtCore.QRect(110, 190, 101, 21))
        self.lineEdit_contourlinspace.setObjectName(_fromUtf8("lineEdit_contourlinspace"))
        self.eigenvaluesButton = QtGui.QPushButton(self.groupBox_Hsm)
        self.eigenvaluesButton.setGeometry(QtCore.QRect(110, 230, 114, 32))
        self.eigenvaluesButton.setAutoDefault(False)
        self.eigenvaluesButton.setObjectName(_fromUtf8("eigenvaluesButton"))
        self.checkBox_withcontour = QtGui.QCheckBox(self.groupBox_Hsm)
        self.checkBox_withcontour.setEnabled(False)
        self.checkBox_withcontour.setGeometry(QtCore.QRect(10, 160, 161, 20))
        self.checkBox_withcontour.setChecked(False)
        self.checkBox_withcontour.setObjectName(_fromUtf8("checkBox_withcontour"))
        self.groupBox_ordering = QtGui.QGroupBox(EnergyDomainDialog)
        self.groupBox_ordering.setEnabled(False)
        self.groupBox_ordering.setGeometry(QtCore.QRect(10, 180, 191, 121))
        self.groupBox_ordering.setObjectName(_fromUtf8("groupBox_ordering"))
        self.orderingButton = QtGui.QPushButton(self.groupBox_ordering)
        self.orderingButton.setGeometry(QtCore.QRect(100, 90, 91, 32))
        self.orderingButton.setAutoDefault(False)
        self.orderingButton.setObjectName(_fromUtf8("orderingButton"))
        self.comboBox_ordering = QtGui.QComboBox(self.groupBox_ordering)
        self.comboBox_ordering.setGeometry(QtCore.QRect(0, 20, 141, 26))
        self.comboBox_ordering.setObjectName(_fromUtf8("comboBox_ordering"))
        self.comboBox_ordering.addItem(_fromUtf8(""))
        self.comboBox_ordering.addItem(_fromUtf8(""))
        self.comboBox_ordering.addItem(_fromUtf8(""))
        self.comboBox_ordering.addItem(_fromUtf8(""))
        self.lineEdit_ordering = QtGui.QLineEdit(self.groupBox_ordering)
        self.lineEdit_ordering.setEnabled(False)
        self.lineEdit_ordering.setGeometry(QtCore.QRect(130, 60, 50, 21))
        self.lineEdit_ordering.setObjectName(_fromUtf8("lineEdit_ordering"))
        self.label_8 = QtGui.QLabel(self.groupBox_ordering)
        self.label_8.setGeometry(QtCore.QRect(70, 60, 51, 16))
        self.label_8.setObjectName(_fromUtf8("label_8"))

        self.retranslateUi(EnergyDomainDialog)
        QtCore.QObject.connect(self.CloseButton, QtCore.SIGNAL(_fromUtf8("clicked()")), EnergyDomainDialog.reject)
        QtCore.QObject.connect(self.lineEdit_QuantumNumber, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.drawButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_hsmprange, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.drawButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_hsmqrange, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.drawButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_hsmgrid, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.drawButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_dim, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.GetEigenButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_pdomain, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.GetEigenButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_qdomain, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.GetEigenButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_ordering, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.orderingButton.animateClick)
        QtCore.QObject.connect(self.lineEdit_contourlinspace, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.drawButton.animateClick)
        QtCore.QMetaObject.connectSlotsByName(EnergyDomainDialog)
        EnergyDomainDialog.setTabOrder(self.comboBox_SelectAlgorism, self.lineEdit_dim)
        EnergyDomainDialog.setTabOrder(self.lineEdit_dim, self.lineEdit_qdomain)
        EnergyDomainDialog.setTabOrder(self.lineEdit_qdomain, self.lineEdit_pdomain)
        EnergyDomainDialog.setTabOrder(self.lineEdit_pdomain, self.comboBox_ordering)
        EnergyDomainDialog.setTabOrder(self.comboBox_ordering, self.lineEdit_ordering)
        EnergyDomainDialog.setTabOrder(self.lineEdit_ordering, self.GetEigenButton)
        EnergyDomainDialog.setTabOrder(self.GetEigenButton, self.lineEdit_hsmgrid)
        EnergyDomainDialog.setTabOrder(self.lineEdit_hsmgrid, self.lineEdit_hsmqrange)
        EnergyDomainDialog.setTabOrder(self.lineEdit_hsmqrange, self.lineEdit_hsmprange)
        EnergyDomainDialog.setTabOrder(self.lineEdit_hsmprange, self.lineEdit_QuantumNumber)
        EnergyDomainDialog.setTabOrder(self.lineEdit_QuantumNumber, self.drawButton)
        EnergyDomainDialog.setTabOrder(self.drawButton, self.CloseButton)

    def retranslateUi(self, EnergyDomainDialog):
        EnergyDomainDialog.setWindowTitle(_translate("EnergyDomainDialog", "Dialog", None))
        self.CloseButton.setText(_translate("EnergyDomainDialog", "Close", None))
        self.groupBox_System.setTitle(_translate("EnergyDomainDialog", "Quantum system setting", None))
        self.GetEigenButton.setText(_translate("EnergyDomainDialog", "get eigen", None))
        self.lineEdit_dim.setText(_translate("EnergyDomainDialog", "50", None))
        self.label.setText(_translate("EnergyDomainDialog", "Hil. dim.", None))
        self.label_3.setText(_translate("EnergyDomainDialog", "p:", None))
        self.label_2.setText(_translate("EnergyDomainDialog", "q:", None))
        self.comboBox_SelectAlgorism.setItemText(0, _translate("EnergyDomainDialog", "Quantum Mapping (SplitOperator)", None))
        self.comboBox_SelectAlgorism.setItemText(1, _translate("EnergyDomainDialog", "Hamiltonian (PositionBase)", None))
        self.comboBox_SelectAlgorism.setItemText(2, _translate("EnergyDomainDialog", "Hamiltonian (HarmonicBase)", None))
        self.groupBox_Hsm.setTitle(_translate("EnergyDomainDialog", "figure setting", None))
        self.label_4.setText(_translate("EnergyDomainDialog", "vq:", None))
        self.label_6.setText(_translate("EnergyDomainDialog", "vp:", None))
        self.label_7.setText(_translate("EnergyDomainDialog", "grid (col,row)", None))
        self.lineEdit_hsmgrid.setText(_translate("EnergyDomainDialog", "50,50", None))
        self.drawButton.setText(_translate("EnergyDomainDialog", "Draw ", None))
        self.label_5.setText(_translate("EnergyDomainDialog", "n:", None))
        self.lineEdit_QuantumNumber.setText(_translate("EnergyDomainDialog", "0", None))
        self.label_9.setText(_translate("EnergyDomainDialog", "min, max, num", None))
        self.lineEdit_contourlinspace.setText(_translate("EnergyDomainDialog", "0,1,10", None))
        self.eigenvaluesButton.setText(_translate("EnergyDomainDialog", "eigen values", None))
        self.checkBox_withcontour.setText(_translate("EnergyDomainDialog", "with energy contour", None))
        self.groupBox_ordering.setTitle(_translate("EnergyDomainDialog", "ordering", None))
        self.orderingButton.setText(_translate("EnergyDomainDialog", "ordering", None))
        self.comboBox_ordering.setItemText(0, _translate("EnergyDomainDialog", "ref. integ. Ham.", None))
        self.comboBox_ordering.setItemText(1, _translate("EnergyDomainDialog", "energy levels", None))
        self.comboBox_ordering.setItemText(2, _translate("EnergyDomainDialog", "q-vari", None))
        self.comboBox_ordering.setItemText(3, _translate("EnergyDomainDialog", "p-vari", None))
        self.lineEdit_ordering.setText(_translate("EnergyDomainDialog", "0", None))
        self.label_8.setText(_translate("EnergyDomainDialog", "center:", None))

