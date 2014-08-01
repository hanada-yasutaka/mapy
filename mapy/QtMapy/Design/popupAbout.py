# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'popupAbout.ui'
#
# Created: Fri Jul 25 11:53:33 2014
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

class Ui_AboutDialog(object):
    def setupUi(self, AboutDialog):
        AboutDialog.setObjectName(_fromUtf8("AboutDialog"))
        AboutDialog.resize(400, 300)
        self.label = QtGui.QLabel(AboutDialog)
        self.label.setGeometry(QtCore.QRect(50, 20, 291, 21))
        self.label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.label.setObjectName(_fromUtf8("label"))
        self.pushButton = QtGui.QPushButton(AboutDialog)
        self.pushButton.setGeometry(QtCore.QRect(260, 250, 114, 32))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.textBrowser = QtGui.QTextBrowser(AboutDialog)
        self.textBrowser.setGeometry(QtCore.QRect(70, 50, 256, 192))
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))

        self.retranslateUi(AboutDialog)
        QtCore.QObject.connect(self.pushButton, QtCore.SIGNAL(_fromUtf8("clicked()")), AboutDialog.close)
        QtCore.QMetaObject.connectSlotsByName(AboutDialog)

    def retranslateUi(self, AboutDialog):
        AboutDialog.setWindowTitle(_translate("AboutDialog", "Dialog", None))
        self.label.setText(_translate("AboutDialog", "QtMapy", None))
        self.pushButton.setText(_translate("AboutDialog", "OK", None))
        self.textBrowser.setHtml(_translate("AboutDialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'.Lucida Grande UI\'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>", None))

