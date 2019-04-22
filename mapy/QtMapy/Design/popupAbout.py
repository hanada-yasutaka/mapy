# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'popupAbout.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_AboutDialog(object):
    def setupUi(self, AboutDialog):
        AboutDialog.setObjectName("AboutDialog")
        AboutDialog.resize(400, 300)
        self.label = QtWidgets.QLabel(AboutDialog)
        self.label.setGeometry(QtCore.QRect(50, 20, 291, 21))
        self.label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.label.setObjectName("label")
        self.pushButton = QtWidgets.QPushButton(AboutDialog)
        self.pushButton.setGeometry(QtCore.QRect(260, 250, 114, 32))
        self.pushButton.setObjectName("pushButton")
        self.textBrowser = QtWidgets.QTextBrowser(AboutDialog)
        self.textBrowser.setGeometry(QtCore.QRect(70, 50, 256, 192))
        self.textBrowser.setObjectName("textBrowser")

        self.retranslateUi(AboutDialog)
        self.pushButton.clicked.connect(AboutDialog.close)
        QtCore.QMetaObject.connectSlotsByName(AboutDialog)

    def retranslateUi(self, AboutDialog):
        _translate = QtCore.QCoreApplication.translate
        AboutDialog.setWindowTitle(_translate("AboutDialog", "Dialog"))
        self.label.setText(_translate("AboutDialog", "QtMapy"))
        self.pushButton.setText(_translate("AboutDialog", "OK"))
        self.textBrowser.setHtml(_translate("AboutDialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'.Lucida Grande UI\'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))

