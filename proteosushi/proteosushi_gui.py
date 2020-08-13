#! /usr/bin/env python3

"""proteosushi_gui.py: displays the GUI and connects it to the main program"""

import os
import sys
from PyQt5 import QtCore
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QPushButton, 
                             QMessageBox, QLineEdit, QLabel, QGroupBox, 
                             QGridLayout, QVBoxLayout, QFileDialog, QCheckBox,
                             QRadioButton, QButtonGroup)
from PyQt5.QtCore import pyqtSlot, QSize
from PyQt5.QtGui import QIcon

from combine_intensities import parse_output, rollup
from parse_proteome import parse_proteome
from proteoSushi_constants import cleave_rules

class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = "ProteoSushi"
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480
        self.centralWidget = QWidget()
        self.setCentralWidget(self.centralWidget)
        self.species_id_dict = self.__make_species_dict()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)


        self.fileChooserLabel = QLabel("Types of Files", self)

        self.maxquant_RB = QRadioButton("MaxQuant", self)
        self.maxquant_RB.toggled.connect(self.check_MQ_RB)
        self.maxquant_RB.setToolTip("Use the txt output folder from the MaxQuant Search Engine")
        self.MQbutton = QPushButton("MaxQuant Output Folder")
        self.MQbutton.setHidden(True)
        self.MQbutton.clicked.connect(self.onClickMQbutton)
        self.maxquantFilepath = QLabel("[Filepath]", self)
        self.maxquantFilepath.setHidden(True)

        self.mascot_RB = QRadioButton("Mascot", self)
        self.mascot_RB.toggled.connect(self.check_mascot_RB)
        self.mascot_RB.setToolTip("Use the output from the Mascot Search Engine")
        self.mascotButton = QPushButton("Mascot Output")
        self.mascotButton.setHidden(True)
        self.mascotButton.clicked.connect(self.onClickMascotButton)
        self.mascotFilepath = QLabel("[Filepath]", self)
        self.mascotFilepath.setHidden(True)

        self.generic_RB = QRadioButton("Generic", self)
        self.generic_RB.toggled.connect(self.check_generic_RB)
        self.generic_RB.setToolTip("Use the output from any other search engine")
        self.genericButton = QPushButton("Generic Output")
        self.genericButton.setHidden(True)
        self.genericButton.clicked.connect(self.onClickGenericButton)
        self.genericFilepath = QLabel("[Filepath]", self)
        self.genericFilepath.setHidden(True)

        self.search_engine_group = QButtonGroup()
        self.search_engine_group.addButton(self.maxquant_RB)
        self.search_engine_group.addButton(self.mascot_RB)
        self.search_engine_group.addButton(self.generic_RB)

        self.PTM_CBs = list()

        self.proteome_filepath_button = QPushButton("Uniprot Proteome FASTA")
        self.proteome_filepath_button.clicked.connect(self.on_click_proteome_button)
        self.proteome_filepath = QLabel("[Filepath]", self)

        self.options_label = QLabel("Options", self)

        self.targetCB = QCheckBox("Use Target Genes", self)
        self.targetCB.stateChanged.connect(self.checkTargetBox)
        self.targetCB.setToolTip("Use a list of genes that will be prioritized given multiple matches")
        self.targetButton = QPushButton("Target Gene List")
        self.targetButton.setHidden(True)
        self.targetButton.clicked.connect(self.onClickTargetButton)
        self.targetFilepath = QLabel("[Filepath]", self)
        self.targetFilepath.setHidden(True)

        self.quant_CB = QCheckBox("Use Quantitation Values", self)
        self.quant_CB.stateChanged.connect(self.check_quant_CB)
        self.quant_CB.setToolTip("Whether to use the Intensity/Quantitation values from the PSMs")
        self.sum_RB = QRadioButton("Sum Peaks", self)
        self.sum_RB.setHidden(True)
        self.average_RB = QRadioButton("Average Peaks", self)
        self.average_RB.setHidden(True)

        self.quant_method_group = QButtonGroup()
        self.quant_method_group.addButton(self.sum_RB)
        self.quant_method_group.addButton(self.average_RB)

        self.uniprot_annot_CB = QCheckBox("Annotate with Uniprot", self)
        self.uniprot_annot_CB.setToolTip("Whether to add in annotation from Uniprot, like subcellular location, secondary structure, etc.")

        self.species_id_label = QLabel("Species ID", self)
        self.species_id_label.setToolTip("The species ID for the species from the data provided (e.g. 9606)")
        self.species_id_edit = QLineEdit(self)
        self.species_id_edit.editingFinished.connect(self.update_species_name)
        self.species_id_edit.setToolTip("The species ID for the species from the data provided (e.g. 9606)")
        self.species_name_label = QLabel("", self)

        self.max_missed_label = QLabel("Max Missed Cleavages", self)
        self.max_missed_label.setToolTip("The maximum allowed missed cleavages for a given peptide")
        self.maxMissedEdit = QLineEdit(self)
        self.maxMissedEdit.setToolTip("The maximum allowed missed cleavages for a given peptide")

        self.fdr_label = QLabel("FDR Threshold", self)
        self.fdr_label.setToolTip("The threshold for pep_expect column for Mascot or PEP column for Maxquant.\nMust be between 0 and 1")
        self.fdrEdit = QLineEdit(self)  # TODO: Error check for a number
        self.fdrEdit.setToolTip("The threshold for pep_expect column for Mascot or PEP column for Maxquant.\nMust be between 0 and 1")

        self.protease_label = QLabel("Protease used in sample digestion", self)
        self.protease_label.setToolTip("The protease used to digest the sample\nExamples include: trypsin/p, trypsin!p, lys-c, asp-n, asp-nc, lys-n")
        self.protease_edit = QLineEdit(self)
        self.protease_edit.setToolTip("The protease used to digest the sample\nExamples include: trypsin/p, trypsin!p, lys-c, asp-n, asp-nc, lys-n")

        self.runButton = QPushButton("Rollup!")
        self.runButton.setHidden(False)
        self.runButton.clicked.connect(self.onClickRunButton)

        self.createGridLayout()
        windowLayout = QVBoxLayout(self.centralWidget)
        windowLayout.addWidget(self.horizontalGroupBox)
        self.setLayout(windowLayout)
        
        self.show()
    
    def check_MQ_RB(self, state):
        if self.maxquant_RB.isChecked():
            self.maxquantFilepath.setHidden(False)
            self.MQbutton.setHidden(False)
        else:
            self.maxquantFilepath.setHidden(True)
            self.MQbutton.setHidden(True)

    def check_mascot_RB(self, state):
        if self.mascot_RB.isChecked():
            self.mascotFilepath.setHidden(False)
            self.mascotButton.setHidden(False)
        else:
            self.mascotFilepath.setHidden(True)
            self.mascotButton.setHidden(True)
    
    def check_generic_RB(self, state):
        if self.generic_RB.isChecked():
            self.genericFilepath.setHidden(False)
            self.genericButton.setHidden(False)
        else:
            self.genericFilepath.setHidden(True)
            self.genericButton.setHidden(True)

    def checkTargetBox(self, state):
        if state == QtCore.Qt.Checked:
            self.targetFilepath.setHidden(False)
            self.targetButton.setHidden(False)
        else:
            self.targetFilepath.setHidden(True)
            self.targetButton.setHidden(True)
    
    def check_quant_CB(self, state):
        if state == QtCore.Qt.Checked:
            self.average_RB.setHidden(False)
            self.sum_RB.setHidden(False)
        else:
            self.average_RB.setHidden(True)
            self.sum_RB.setHidden(True)
    
    def update_species_name(self):
        if not self.species_id_edit.text(): 
            self.species_name_label.setStyleSheet("background-color : white")
            self.species_name_label.setText("")
        elif self.species_id_edit.text() in self.species_id_dict:
            self.species_name_label.setStyleSheet("background-color : white")
            self.species_name_label.setText(self.species_id_dict[self.species_id_edit.text()][2:])
        else:
            self.species_name_label.setText("Not a valid species ID")
            self.species_name_label.setStyleSheet("background-color : red")

    def checkRunButton(self, MQcheckState, MQfilepath: str, mascotCheckState, mascotFilepath: str, genericCheckState, genericFilepath: str):
        """Checks to see whether the run button is available. must have a search engine output and missed cleavages, etc
        Arguments:
            MQcheckState {QtCore.Qt.Checked} -- whether the box is checked
            MQfilepath {str} -- the filepath for maxquant output
            mascotCheckState {QtCore.Qt.Checked} -- whether the box is checked
            mascotFilepath {str} -- the filepath for maxquant output
            genericCheckState {QtCore.Qt.Checked} -- whether the box is checked
            genericFilepath {str} -- the filepath for maxquant output
        """
        # TODO: update to require species, etc.
        if (((MQcheckState == QtCore.Qt.Checked and MQfilepath != "") or
            (mascotCheckState == QtCore.Qt.Checked and mascotFilepath != "") or
            (genericCheckState == QtCore.Qt.Checked and genericFilepath != "")) and
            os.path.exists(self.proteome_filepath.text())):
            self.runButton.setHidden(False)
        else:
            self.runButton.setHidden(True)


    def __make_species_dict(self):
        species_id_dict = dict()
        with open("spec_list_fixed.tsv", 'r') as spec_list:
            for line in spec_list:
                split_line = line.strip().split('\t')
                species_id_dict[split_line[0]] = split_line[1]
        return species_id_dict



    def openCSVFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","CSV Files (*.csv);;All Files (*)", options=options)
        return fileName

    def openTXTFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","TXT Files (*.txt);;All Files (*)", options=options)
        return fileName

    def open_directory_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        #dir_name, _ = QFileDialog.getExistingDirectory(self,"QFileDialog.getExistingDirectory()", "","TXT Files (*.txt);;All Files (*)", options=options)
        dir_name = str(QFileDialog.getExistingDirectory(self, "Select Directory", options=options))
        return dir_name
    
    def openFASTAFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","FASTA Files (*.fasta);;All Files (*)", options=options)
        return fileName


    def convert_PTM_list(self, PTM_CBs:list) -> list:
        """changes the list of checkboxes into chosen ptms for analysis

        Arguments:
            PTM_CBs {list} -- a list of QCheckBoxes for each possible PTM
        Returns:
            list -- a list of strings of PTM names for analysis
        """
        PTM_str_list = list()
        for ptm_cb in PTM_CBs:
            if ptm_cb.isChecked():
                PTM_str_list.append(ptm_cb.text())
        return PTM_str_list

    @pyqtSlot()
    def onClickMQbutton(self):
        self.statusBar().showMessage("Choose the MaxQuant output FOLDER")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.open_directory_dialog()
        if self.maxquantFilepath.text() != "":
            self.maxquantFilepath.setText(filename)
            missed_cleavages, enzyme, PTMs = parse_output("maxquant", filename)
            self.maxMissedEdit.setText(str(missed_cleavages))
            self.protease_edit.setText(enzyme)
            #Remove the previous PTM checkboxes (if there were any)
            if self.PTM_CBs != []:
                for cb in self.PTM_CBs:
                    self.layout.removeWidget(cb)
                    cb.deleteLater()

            # Adds the new checkboxes to the list
            self.PTM_CBs = []
            for ptm in PTMs:
                self.PTM_CBs.append(QCheckBox(ptm, self))
            # Removes the label
            i = 0
            if not self.layout.itemAtPosition(5, i) is None:
                label_to_remove = self.layout.itemAtPosition(5, i).widget()
                self.layout.removeWidget(label_to_remove)
                label_to_remove.deleteLater()

            # Inserts the new PTM checkboxes
            self.layout.addWidget(QLabel("PTMs for Analysis"), 5, i)
            for widget in self.PTM_CBs:
                self.layout.addWidget(widget, 6, i)
                i += 1
        self.statusBar().showMessage("")
            
        '''
        buttonReply = QMessageBox.warning(self, "ERROR", "Error detected", QMessageBox.Ok, QMessageBox.Ok)
        if buttonReply == QMessageBox.Yes:
            print('Yes clicked.')
        else:
            print('No clicked.')
        '''
    
    @pyqtSlot()
    def onClickMascotButton(self):
        self.statusBar().showMessage("Choose the Mascot output file")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openCSVFileNameDialog()
        if self.mascotFilepath.text() != "":
            self.mascotFilepath.setText(filename)
            missed_cleavages, enzyme, PTMs = parse_output("mascot", filename)
            self.maxMissedEdit.setText(str(missed_cleavages))
            self.protease_edit.setText(enzyme)
            #Remove the previous PTM checkboxes (if there were any)
            if self.PTM_CBs != []:
                for cb in self.PTM_CBs:
                    self.layout.removeWidget(cb)
                    cb.deleteLater()

            # Adds the new checkboxes to the list
            self.PTM_CBs = []
            for ptm in PTMs:
                self.PTM_CBs.append(QCheckBox(ptm, self))
            # Removes the label
            i = 0
            if not self.layout.itemAtPosition(5, i) is None:
                label_to_remove = self.layout.itemAtPosition(5, i).widget()
                self.layout.removeWidget(label_to_remove)
                label_to_remove.deleteLater()
            #self.horizontalGroupBox.setLayout(self.layout)
            
            # Inserts the new PTM checkboxes
            self.layout.addWidget(QLabel("PTMs for Analysis"), 5, i)
            for widget in self.PTM_CBs:
                #if not self.layout.itemAtPosition(6, i) is None:
                #    self.layout.removeWidget(self.layout.itemAtPosition(6, i).widget())
                self.layout.addWidget(widget, 6, i)
                i += 1
            #self.horizontalGroupBox.setLayout(self.layout)
        self.statusBar().showMessage("")
    
    @pyqtSlot()
    def onClickGenericButton(self):
        self.statusBar().showMessage("Choose the Search Engine output")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openCSVFileNameDialog()
        if self.genericFilepath.text() != "":
            self.genericFilepath.setText(filename)
            missed_cleavages, enzyme, PTMs = parse_output("generic", filename)
            #Remove the previous PTM checkboxes (if there were any)
            if self.PTM_CBs != []:
                for cb in self.PTM_CBs:
                    self.layout.removeWidget(cb)
                    cb.deleteLater()

            # Adds the new checkboxes to the list
            self.PTM_CBs = []
            for ptm in PTMs:
                self.PTM_CBs.append(QCheckBox(ptm, self))
            # Removes the label
            i = 0
            if not self.layout.itemAtPosition(5, i) is None:
                label_to_remove = self.layout.itemAtPosition(5, i).widget()
                self.layout.removeWidget(label_to_remove)
                label_to_remove.deleteLater()

            # Inserts the new PTM checkboxes
            self.layout.addWidget(QLabel("PTMs for Analysis"), 5, i)
            for widget in self.PTM_CBs:
                self.layout.addWidget(widget, 6, i)
                i += 1
        self.statusBar().showMessage("")

    @pyqtSlot()
    def on_click_proteome_button(self):
        self.statusBar().showMessage("Choose the Uniprot Proteome FASTA file")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openFASTAFileNameDialog()
        if self.proteome_filepath.text() != "":
            # Parse the file to check it is valid and grab the species
            was_error, species_ID = parse_proteome(filename)
            if was_error:
                self.statusBar().showMessage("ERROR: Invalid proteome FASTA file!")
                self.statusBar().setStyleSheet("background-color : red")
            else:
                self.proteome_filepath.setText(filename)
                if not species_ID is None:
                    self.species_id_edit.setText(species_ID)
        self.statusBar().showMessage("")

    @pyqtSlot()
    def onClickTargetButton(self):
        self.statusBar().showMessage("Choose the Target Gene file")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openTXTFileNameDialog()
        if self.targetFilepath.text() != "":
            self.targetFilepath.setText(filename)
        self.statusBar().showMessage("")

    @pyqtSlot()
    def onClickRunButton(self):
        if ((self.maxquant_RB.isChecked() and os.path.exists(self.maxquantFilepath.text())) or
            (self.mascot_RB.isChecked() and os.path.exists(self.mascotFilepath.text())) or
            (self.generic_RB.isChecked() and os.path.exists(self.genericFilepath.text()))): 
            if not os.path.exists(self.proteome_filepath.text()):
                self.statusBar().showMessage("ERROR: Missing Proteome FASTA file!")
                self.statusBar().setStyleSheet("background-color : red")
                return

            try:  # Check to see if the provided number of allowed missed cleavages is legal
                max_missed = int(self.maxMissedEdit.text())
                if max_missed < 0:
                    self.statusBar().showMessage("ERROR: Number of max allowed missed cleavages cannot be less than 0!")
                    self.statusBar().setStyleSheet("background-color : red")
                    return
            except ValueError:
                self.statusBar().showMessage("ERROR: Invalid number for max allowed missed cleavages!")
                self.statusBar().setStyleSheet("background-color : red")
                return
            
            try:  # Checks to see if the provided FDR threshold is legal
                if self.fdrEdit.text() != "":
                    fdr = float(self.fdrEdit.text())
                    if fdr < 0 or fdr > 1:
                        self.statusBar().showMessage("ERROR: FDR must be between 0 and 1!")
                        self.statusBar().setStyleSheet("background-color : red")
                        return
                else:
                    fdr = None
            except ValueError:
                self.statusBar().showMessage("ERROR: FDR must be a number")
                self.statusBar().setStyleSheet("background-color : red")
                return
            # If quantitation will be used, but the user has not selected a method to combine
            if self.quant_CB.isChecked() and not self.sum_RB.isChecked() and not self.average_RB.isChecked():
                self.statusBar().showMessage("ERROR: Must choose either Sum or Average for the quant method!")
                self.statusBar().setStyleSheet("background-color : red")
                return
            elif self.quant_CB.isChecked() and self.sum_RB.isChecked():
                combine_method = "sum"
            elif self.quant_CB.isChecked() and self.average_RB.isChecked():
                combine_method = "average"
            else:
                combine_method = ""
            
            # If a valid species ID has been entered in
            species_id = self.species_id_edit.text()
            #print(species_id)
            if len(species_id) > 0 and not species_id in self.species_id_dict:
                self.statusBar().showMessage("ERROR: Not a valid species ID!")
                self.statusBar().setStyleSheet("background-color : red")
                return

                
            if self.protease_edit.text() in cleave_rules:
                #self.statusBar().showMessage("Analysis Running!")  # TODO: consider running on a different thread
                
                # If the maxquant option was chosen, it sends that info to be run
                if self.maxquant_RB.isChecked() and os.path.exists(self.maxquantFilepath.text()):
                    rollup("maxquant", 
                                    self.maxquantFilepath.text(), 
                                    self.targetCB.isChecked(),  # Whether target will be used
                                    self.targetFilepath.text(),
                                    int(self.maxMissedEdit.text()),
                                    self.protease_edit.text(),
                                    fdr,
                                    self.quant_CB.isChecked(),
                                    self.convert_PTM_list(self.PTM_CBs),
                                    self.proteome_filepath.text(),
                                    combine_method,
                                    self.uniprot_annot_CB.isChecked(),
                                    species_id)
                    self.statusBar().showMessage("Analysis Complete!")
                    self.statusBar().setStyleSheet("background-color : green")
                elif self.mascot_RB.isChecked() and os.path.exists(self.mascotFilepath.text()):
                    rollup("mascot", 
                                    self.mascotFilepath.text(), 
                                    self.targetCB.isChecked(),  # Whether target will be used
                                    self.targetFilepath.text(),
                                    int(self.maxMissedEdit.text()),
                                    self.protease_edit.text(),
                                    fdr,
                                    self.quant_CB.isChecked(),
                                    self.convert_PTM_list(self.PTM_CBs),
                                    self.proteome_filepath.text(),
                                    combine_method,
                                    self.uniprot_annot_CB.isChecked(),
                                    species_id)
                    self.statusBar().showMessage("Analysis Complete!")
                    self.statusBar().setStyleSheet("background-color : green")
                elif self.generic_RB.isChecked() and os.path.exists(self.genericFilepath.text()):
                    rollup("generic", 
                                    self.genericFilepath.text(), 
                                    self.targetCB.isChecked(),  # Whether target will be used
                                    self.targetFilepath.text(),
                                    int(self.maxMissedEdit.text()),
                                    self.protease_edit.text(),
                                    fdr,
                                    self.quant_CB.isChecked(),
                                    self.convert_PTM_list(self.PTM_CBs),
                                    self.proteome_filepath.text(),
                                    combine_method,
                                    self.uniprot_annot_CB.isChecked(),
                                    species_id)
                    self.statusBar().showMessage("Analysis Complete!")
                    self.statusBar().setStyleSheet("background-color : green")
                else:
                    self.statusBar().showMessage("ERROR: Missing Search Engine Output!")
                    self.statusBar().setStyleSheet("background-color : red")
            else:
                self.statusBar().showMessage("ERROR: Protease provided is not supported!")
                self.statusBar().setStyleSheet("background-color : red")
        else:
            self.statusBar().showMessage("ERROR: Missing Search Engine Output!")
            self.statusBar().setStyleSheet("background-color : red")




    def createGridLayout(self):
        self.horizontalGroupBox = QGroupBox("")
        self.layout = QGridLayout()
        self.layout.setRowStretch(0, 1)
        self.layout.setRowStretch(9, 1)

        row = 1
        self.layout.addWidget(self.fileChooserLabel, row, 0)
        row += 1
        self.layout.addWidget(self.maxquant_RB, row, 0)
        self.layout.addWidget(self.MQbutton, row, 1)
        self.layout.addWidget(self.maxquantFilepath, row, 2)
        row += 1
        self.layout.addWidget(self.mascot_RB, row, 0)
        self.layout.addWidget(self.mascotButton, row, 1)
        self.layout.addWidget(self.mascotFilepath, row, 2)
        row += 1
        self.layout.addWidget(self.generic_RB, row, 0)
        self.layout.addWidget(self.genericButton, row, 1)
        self.layout.addWidget(self.genericFilepath, row, 2)
        row += 1
        # This is where the PTM label goes
        row += 1
        # This is where the PTM checkboxes go
        row += 1
        self.layout.addWidget(QLabel("Proteome File", self), row, 0)
        self.layout.addWidget(self.proteome_filepath_button, row, 1)
        self.layout.addWidget(self.proteome_filepath, row, 2)
        row += 1
        self.layout.addWidget(self.options_label, row, 0)
        row += 1
        self.layout.addWidget(self.targetCB, row, 0)
        self.layout.addWidget(self.targetButton, row, 1)
        self.layout.addWidget(self.targetFilepath, row, 2)
        row += 1
        self.layout.addWidget(self.quant_CB, row, 0)
        self.layout.addWidget(self.sum_RB, row, 1)
        self.layout.addWidget(self.average_RB, row, 2)
        row += 1
        self.layout.addWidget(self.uniprot_annot_CB, row, 0)
        row += 1
        self.layout.addWidget(self.species_id_label, row, 0)
        self.layout.addWidget(self.species_id_edit, row, 1)
        self.layout.addWidget(self.species_name_label, row, 2)
        row += 1
        self.layout.addWidget(self.max_missed_label, row, 0)
        self.layout.addWidget(self.maxMissedEdit, row, 1)
        row += 1
        self.layout.addWidget(self.fdr_label, row, 0)
        self.layout.addWidget(self.fdrEdit, row, 1)
        row += 1
        self.layout.addWidget(self.protease_label, row, 0)
        self.layout.addWidget(self.protease_edit, row, 1)
        row += 1
        self.layout.addWidget(self.runButton, row, 0)

        self.horizontalGroupBox.setLayout(self.layout)


def run_gui():
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())

if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
