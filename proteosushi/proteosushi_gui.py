#! /usr/bin/env python3

"""proteosushi_gui.py: displays the GUI and connects it to the main program"""

import importlib.resources as pkg_resources
import os
import sys
from PyQt5 import QtCore
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QPushButton, 
                             QMessageBox, QLineEdit, QLabel, QGroupBox, 
                             QGridLayout, QVBoxLayout, QFileDialog, QCheckBox,
                             QRadioButton, QButtonGroup, QComboBox)
from PyQt5.QtCore import pyqtSlot, QSize, Qt, QRunnable, QObject, QThreadPool, pyqtSignal
from PyQt5.QtGui import QIcon

from .combine_intensities import parse_output, rollup
from . import lib
from .parse_proteome import parse_proteome
from .proteoSushi_constants import cleave_rules


class WorkerSignals(QObject):
    """Class that contains the signal sent at the end of the worker thread"""
    end = pyqtSignal()


class Worker(QRunnable):
    """Class for the backend worker thread"""

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()    


    @pyqtSlot()
    def run(self):
        """Initializes the runner function"""
        # Retrieve args/kwargs here; and fire processing using them
        self.fn(*self.args, **self.kwargs)
        self.signals.end.emit()

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
        self.setWindowTitle("ProteoSushi")
        self.setWindowIcon(QIcon("ProteoSushi_icon.png"))
        self.species_id_dict = self.__make_species_dict()
        self.species_name_dict = {v:k for k, v in self.species_id_dict.items()}
        self.threadpool = QThreadPool()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)


        self.file_chooser_label = QLabel("Types of ProteoSushi Files", self)

        self.maxquant_RB = QRadioButton("MaxQuant", self)
        self.maxquant_RB.toggled.connect(self.check_MQ_RB)
        self.maxquant_RB.setToolTip("Use the txt output folder from the MaxQuant Search Engine")
        self.maxquant_button = QPushButton("MaxQuant txt Folder")
        self.maxquant_button.setHidden(True)
        self.maxquant_button.clicked.connect(self.onClickMQbutton)
        self.maxquant_filepath = QLabel("[Filepath]", self)
        self.maxquant_filepath.setHidden(True)

        self.mascot_RB = QRadioButton("Mascot", self)
        self.mascot_RB.toggled.connect(self.check_mascot_RB)
        self.mascot_RB.setToolTip("Use the output from the Mascot Search Engine")
        self.mascot_button = QPushButton("Mascot File")
        self.mascot_button.setHidden(True)
        self.mascot_button.clicked.connect(self.onClickMascotButton)
        self.mascot_filepath = QLabel("[Filepath]", self)
        self.mascot_filepath.setHidden(True)

        self.generic_RB = QRadioButton("Generic", self)
        self.generic_RB.toggled.connect(self.check_generic_RB)
        self.generic_RB.setToolTip("Use the output from any other search engine")
        self.generic_button = QPushButton("Generic File")
        self.generic_button.setHidden(True)
        self.generic_button.clicked.connect(self.onClickGenericButton)
        self.generic_filepath = QLabel("[Filepath]", self)
        self.generic_filepath.setHidden(True)

        self.search_engine_group = QButtonGroup()
        self.search_engine_group.addButton(self.maxquant_RB)
        self.search_engine_group.addButton(self.mascot_RB)
        self.search_engine_group.addButton(self.generic_RB)

        self.PTM_CBs = list()

        self.proteome_filepath_button = QPushButton("Uniprot Proteome FASTA")
        self.proteome_filepath_button.clicked.connect(self.on_click_proteome_button)
        self.proteome_filepath = QLabel("[Filepath]", self)

        self.output_filepath_button = QPushButton("Output Name and Location")
        self.output_filepath_button.clicked.connect(self.on_click_output_button)
        self.output_filepath = QLabel("[Filepath]", self)

        self.options_label = QLabel("Options", self)

        self.target_checkbox = QCheckBox("Use Target Genes", self)
        self.target_checkbox.stateChanged.connect(self.checkTargetBox)
        self.target_checkbox.setToolTip("Use a list of genes that will be prioritized given multiple matches")
        self.target_button = QPushButton("Target Gene List")
        self.target_button.setHidden(True)
        self.target_button.clicked.connect(self.onClickTargetButton)
        self.target_filepath = QLabel("[Filepath]", self)
        self.target_filepath.setHidden(True)

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
        self.max_missed_edit = QLineEdit(self)
        self.max_missed_edit.setToolTip("The maximum allowed missed cleavages for a given peptide")

        self.fdr_label = QLabel("FDR Threshold", self)
        self.fdr_label.setToolTip("[OPTIONAL] The threshold for pep_expect column for Mascot or PEP column for Maxquant.\nMust be between 0 and 1, but can be left blank.")
        self.fdr_edit = QLineEdit(self)  # TODO: Error check for a number
        self.fdr_edit.setToolTip("[OPTIONAL] The threshold for pep_expect column for Mascot or PEP column for Maxquant.\nMust be between 0 and 1, but can be left blank.")

        self.protease_label = QLabel("Protease used in sample digestion", self)
        self.protease_label.setToolTip("The protease used to digest the sample\nExamples include: trypsin/p, trypsin!p, lys-c, asp-n, asp-nc, lys-n")
        self.protease_combo_box = QComboBox()
        self.protease_combo_box.setToolTip("The protease used to digest the sample\nExamples include: trypsin/p, trypsin!p, lys-c, asp-n, asp-nc, lys-n")
        self.protease_combo_box.addItems(["trypsin/p", "trypsin!p", "lys-c", "asp-n", "asp-nc", "lys-n"])
        self.protease_combo_box.setEditable(True)


        self.run_button = QPushButton("Rollup!")
        self.run_button.setHidden(False)
        self.run_button.clicked.connect(self.onClickRunButton)

        self.createGridLayout()
        windowLayout = QVBoxLayout(self.centralWidget)
        windowLayout.addWidget(self.horizontalGroupBox)
        windowLayout.setAlignment(Qt.AlignTop)
        self.setLayout(windowLayout)
        
        self.show()
    
    def check_MQ_RB(self, state):
        if self.maxquant_RB.isChecked():
            self.maxquant_filepath.setHidden(False)
            self.maxquant_button.setHidden(False)
        else:
            self.maxquant_filepath.setHidden(True)
            self.maxquant_button.setHidden(True)

    def check_mascot_RB(self, state):
        if self.mascot_RB.isChecked():
            self.mascot_filepath.setHidden(False)
            self.mascot_button.setHidden(False)
        else:
            self.mascot_filepath.setHidden(True)
            self.mascot_button.setHidden(True)
    
    def check_generic_RB(self, state):
        if self.generic_RB.isChecked():
            self.generic_filepath.setHidden(False)
            self.generic_button.setHidden(False)
        else:
            self.generic_filepath.setHidden(True)
            self.generic_button.setHidden(True)

    def checkTargetBox(self, state):
        if state == QtCore.Qt.Checked:
            self.target_filepath.setHidden(False)
            self.target_button.setHidden(False)
        else:
            self.target_filepath.setHidden(True)
            self.target_button.setHidden(True)
    
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
        elif f"N={self.species_id_edit.text()}" in self.species_name_dict:
            self.species_name_label.setStyleSheet("background-color : white")
            self.species_name_label.setText(self.species_id_edit.text())
            self.species_id_edit.setText(self.species_name_dict[f"N={self.species_id_edit.text()}"])
        else:
            self.species_name_label.setText("Not a valid species ID")
            self.species_name_label.setStyleSheet("background-color : red")

    def checkRunButton(self, MQcheckState, MQfilepath: str, mascotCheckState, mascot_filepath: str, genericCheckState, generic_filepath: str):
        """Checks to see whether the run button is available. must have a search engine output and missed cleavages, etc
        Arguments:
            MQcheckState {QtCore.Qt.Checked} -- whether the box is checked
            MQfilepath {str} -- the filepath for maxquant output
            mascotCheckState {QtCore.Qt.Checked} -- whether the box is checked
            mascot_filepath {str} -- the filepath for maxquant output
            genericCheckState {QtCore.Qt.Checked} -- whether the box is checked
            generic_filepath {str} -- the filepath for maxquant output
        """
        # TODO: update to require species, etc.
        if (((MQcheckState == QtCore.Qt.Checked and os.path.exists(MQfilepath)) or
            (mascotCheckState == QtCore.Qt.Checked and os.path.exists(mascot_filepath)) or
            (genericCheckState == QtCore.Qt.Checked and os.path.exists(generic_filepath))) and
            os.path.exists(self.proteome_filepath.text())):
            self.run_button.setHidden(False)
        else:
            self.run_button.setHidden(True)


    def __make_species_dict(self):
        species_id_dict = dict()
        with pkg_resources.open_text(lib, "spec_list_fixed.tsv") as spec_list:
        #with pkg_resources.open_text(__package__, "spec_list_fixed.tsv") as spec_list:
            for line in spec_list:
                split_line = line.strip().split('\t')
                species_id_dict[split_line[0]] = split_line[1]
        return species_id_dict



    def openCSVFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Choose the Search Engine Output", os.getcwd(),"CSV Files (*.csv);;All Files (*)", options=options)
        return fileName

    def openTXTFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Choose the Gene Name File", os.getcwd(),"TXT Files (*.txt);;All Files (*)", options=options)
        return fileName

    def open_directory_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dir_name = str(QFileDialog.getExistingDirectory(self, "Select Directory", os.getcwd(), options=options))
        return dir_name
    
    def openFASTAFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Choose the Proteome FASTA file", os.getcwd(),"FASTA Files (*.fasta);;All Files (*)", options=options)
        return fileName
    
    def write_output_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Choose the name and location of the output", os.getcwd(),"CSV Files (*.csv);;All Files (*)", options=options)
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
        #if self.maxquant_filepath.text() != "":
        if os.path.exists(filename):
            self.maxquant_filepath.setText(filename)
            missed_cleavages, enzyme, PTMs = parse_output("maxquant", filename)
            self.max_missed_edit.setText(str(missed_cleavages))
            self.protease_combo_box.setEditText(enzyme)
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
    def onClickMascotButton(self):
        self.statusBar().showMessage("Choose the Mascot output file")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openCSVFileNameDialog()
        #if self.mascot_filepath.text() != "":
        if os.path.exists(filename):
            self.mascot_filepath.setText(filename)
            missed_cleavages, enzyme, PTMs = parse_output("mascot", filename)
            if missed_cleavages == -5:
                self.statusBar().showMessage("Invalid Mascot file")
                self.statusBar().setStyleSheet("background-color : red")
                self.mascot_filepath.setText("")
                return
            self.max_missed_edit.setText(str(missed_cleavages))
            self.protease_combo_box.setEditText(enzyme)
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
    def onClickGenericButton(self):
        self.statusBar().showMessage("Choose the Search Engine output")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openCSVFileNameDialog()
        #if self.generic_filepath.text() != "":
        if os.path.exists(filename):
            self.generic_filepath.setText(filename)
            missed_cleavages, enzyme, PTMs = parse_output("generic", filename)
            if missed_cleavages == -3:
                self.statusBar().showMessage("A sequence in the Peptide Modified Sequence column is missing PTMs")
                self.statusBar().setStyleSheet("background-color : red")
                return
            if missed_cleavages == -4:
                self.statusBar().showMessage("Invalid file")
                self.statusBar().setStyleSheet("background-color : red")
                self.generic_filepath.setText("")
                return
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
        #if self.proteome_filepath.text() != "":
        if os.path.exists(filename):
            # Parse the file to check it is valid and grab the species
            was_error, species_ID = parse_proteome(filename)
            if was_error:
                self.statusBar().showMessage("ERROR: Invalid proteome FASTA file!")
                self.statusBar().setStyleSheet("background-color : red")
            else:
                self.proteome_filepath.setText(filename)
                if not species_ID is None:
                    self.species_id_edit.setText(species_ID)
                    self.update_species_name()
        self.statusBar().showMessage("")

    @pyqtSlot()
    def on_click_output_button(self):
        self.statusBar().showMessage("Choose the ProteoSushi output file location")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.write_output_dialog()
        self.output_filepath.setText(filename)
        self.statusBar().showMessage("")

    @pyqtSlot()
    def onClickTargetButton(self):
        self.statusBar().showMessage("Choose the Target Gene file")
        self.statusBar().setStyleSheet("background-color : white")
        filename = self.openTXTFileNameDialog()
        if self.target_filepath.text() != "":
            self.target_filepath.setText(filename)
        self.statusBar().showMessage("")

    def __run_async(self, fdr, combine_method, species_id):
        """Function that runs the backend on a separate thread"""
        # If the maxquant option was chosen, it sends that info to be run
        if self.maxquant_RB.isChecked() and os.path.exists(self.maxquant_filepath.text()):
            self.statusBar().showMessage("Analysis in Progress")
            self.statusBar().setStyleSheet("background-color : white")
            output = rollup("maxquant", 
                            self.maxquant_filepath.text(), 
                            self.target_checkbox.isChecked(),  # Whether target will be used
                            self.target_filepath.text(),
                            int(self.max_missed_edit.text()),
                            self.protease_combo_box.currentText(),
                            fdr,
                            self.quant_CB.isChecked(),
                            self.convert_PTM_list(self.PTM_CBs),
                            self.proteome_filepath.text(),
                            combine_method,
                            self.uniprot_annot_CB.isChecked(),
                            species_id,
                            self.output_filepath.text())
            # If there is a 502 proxy error (server side error)
            if self.uniprot_annot_CB.isChecked() and output == 502:
                self.statusBar().showMessage("ERROR: Uniprot server error! Please try again later.")
                self.statusBar().setStyleSheet("background-color : red")
                return
            self.statusBar().showMessage("Analysis Complete!")
            self.statusBar().setStyleSheet("background-color : green")
            print("\033[92m {}\033[00m".format("Analysis Complete!"))
            #sys.exit()
        elif self.mascot_RB.isChecked() and os.path.exists(self.mascot_filepath.text()):
            self.statusBar().showMessage("Analysis in Progress")
            self.statusBar().setStyleSheet("background-color : white")
            output = rollup("mascot", 
                            self.mascot_filepath.text(), 
                            self.target_checkbox.isChecked(),  # Whether target will be used
                            self.target_filepath.text(),
                            int(self.max_missed_edit.text()),
                            self.protease_combo_box.currentText(),
                            fdr,
                            self.quant_CB.isChecked(),
                            self.convert_PTM_list(self.PTM_CBs),
                            self.proteome_filepath.text(),
                            combine_method,
                            self.uniprot_annot_CB.isChecked(),
                            species_id,
                            self.output_filepath.text())
            # If the mascot file has no intensity values and user tried to analyze them
            if self.quant_CB.isChecked() and output == 2:
                self.statusBar().showMessage("ERROR: Mascot file has no detectable intensity values!")
                self.statusBar().setStyleSheet("background-color : red")
                return
            # If there is a 502 proxy error (server side error)
            if self.uniprot_annot_CB.isChecked() and output == 502:
                self.statusBar().showMessage("ERROR: Uniprot server error! Please try again later.")
                self.statusBar().setStyleSheet("background-color : red")
                return
            self.statusBar().showMessage("Analysis Complete!")
            self.statusBar().setStyleSheet("background-color : green")
            print("\033[92m {}\033[00m".format("Analysis Complete!"))
            #sys.exit()
        elif self.generic_RB.isChecked() and os.path.exists(self.generic_filepath.text()):
            self.statusBar().showMessage("Analysis in Progress")
            self.statusBar().setStyleSheet("background-color : white")
            output = rollup("generic", 
                            self.generic_filepath.text(), 
                            self.target_checkbox.isChecked(),  # Whether target will be used
                            self.target_filepath.text(),
                            int(self.max_missed_edit.text()),
                            self.protease_combo_box.currentText(),
                            fdr,
                            self.quant_CB.isChecked(),
                            self.convert_PTM_list(self.PTM_CBs),
                            self.proteome_filepath.text(),
                            combine_method,
                            self.uniprot_annot_CB.isChecked(),
                            species_id,
                            self.output_filepath.text())
            # If there is a 502 proxy error (server side error)
            if self.uniprot_annot_CB.isChecked() and output == 502:
                self.statusBar().showMessage("ERROR: Uniprot server error! Please try again later.")
                self.statusBar().setStyleSheet("background-color : red")
                return
            self.statusBar().showMessage("Analysis Complete!")
            self.statusBar().setStyleSheet("background-color : green")
            print("\033[92m {}\033[00m".format("Analysis Complete!"))
            #sys.exit()
        else:
            self.statusBar().showMessage("ERROR: Missing Search Engine Output!")
            self.statusBar().setStyleSheet("background-color : red")
        #sys.exit()

    def __cut_thread(self):
        """Called once the backend thread finishes"""
        # Enables the buttons again in case I get multiple analyses working
        self.run_button.setEnabled(True)
        self.proteome_filepath_button.setEnabled(True)
        self.output_filepath_button.setEnabled(True)
        self.maxquant_button.setEnabled(True)
        self.mascot_button.setEnabled(True)
        self.generic_button.setEnabled(True)
        #sys.exit()
        return

    @pyqtSlot()
    def onClickRunButton(self):
        if ((self.maxquant_RB.isChecked() and os.path.exists(self.maxquant_filepath.text())) or
            (self.mascot_RB.isChecked() and os.path.exists(self.mascot_filepath.text())) or
            (self.generic_RB.isChecked() and os.path.exists(self.generic_filepath.text()))): 
            if not os.path.exists(self.proteome_filepath.text()):
                self.statusBar().showMessage("ERROR: Missing Proteome FASTA file!")
                self.statusBar().setStyleSheet("background-color : red")
                return

            # If the user didn't choose a name or location for the file
            if self.output_filepath.text() == "":
                self.statusBar().showMessage("ERROR: ProteoSushi output name and location not chosen!")
                self.statusBar().setStyleSheet("background-color : red")
                return
            
            if os.path.isdir(self.output_filepath.text()):
                self.output_filepath.setText(os.path.join(self.output_filepath.text(), "proteosushi_output.csv"))

            if self.output_filepath.text() == "[Filepath]":
                self.output_filepath.setText(os.path.join(os.getcwd(), "proteosushi_output.csv"))

            if self.output_filepath.text()[-4:] != ".csv":
                self.output_filepath.setText(self.output_filepath.text() + ".csv")

            try:  # Check to see if the provided number of allowed missed cleavages is legal
                max_missed = int(self.max_missed_edit.text())
                if max_missed < 0:
                    self.statusBar().showMessage("ERROR: Number of max allowed missed cleavages cannot be less than 0!")
                    self.statusBar().setStyleSheet("background-color : red")
                    return
            except ValueError:
                self.statusBar().showMessage("ERROR: Invalid number for max allowed missed cleavages!")
                self.statusBar().setStyleSheet("background-color : red")
                return
            
            try:  # Checks to see if the provided FDR threshold is legal
                if self.fdr_edit.text() != "":
                    fdr = float(self.fdr_edit.text())
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

                
            if self.protease_combo_box.currentText() in cleave_rules:
                # Disable all buttons before running to keep analysis confined
                self.run_button.setEnabled(False)
                self.proteome_filepath_button.setEnabled(False)
                self.output_filepath_button.setEnabled(False)
                self.maxquant_button.setEnabled(False)
                self.mascot_button.setEnabled(False)
                self.generic_button.setEnabled(False)
                # Starts a new thread for the backend analysis
                worker = Worker(self.__run_async, fdr, combine_method, species_id)
                worker.signals.end.connect(self.__cut_thread)
                self.threadpool.start(worker)
                
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
        self.layout.addWidget(self.file_chooser_label, row, 0)
        row += 1
        self.layout.addWidget(self.maxquant_RB, row, 0)
        self.layout.addWidget(self.maxquant_button, row, 1)
        self.layout.addWidget(self.maxquant_filepath, row, 2)
        row += 1
        self.layout.addWidget(self.mascot_RB, row, 0)
        self.layout.addWidget(self.mascot_button, row, 1)
        self.layout.addWidget(self.mascot_filepath, row, 2)
        row += 1
        self.layout.addWidget(self.generic_RB, row, 0)
        self.layout.addWidget(self.generic_button, row, 1)
        self.layout.addWidget(self.generic_filepath, row, 2)
        row += 1
        # This is where the PTM label goes
        row += 1
        # This is where the PTM checkboxes go
        row += 1
        self.layout.addWidget(QLabel("Proteome File", self), row, 0)
        self.layout.addWidget(self.proteome_filepath_button, row, 1)
        self.layout.addWidget(self.proteome_filepath, row, 2)
        row += 1
        self.layout.addWidget(QLabel("Output Name & Location", self), row, 0)
        self.layout.addWidget(self.output_filepath_button, row, 1)
        self.layout.addWidget(self.output_filepath, row, 2)
        row += 1
        self.layout.addWidget(self.options_label, row, 0)
        row += 1
        self.layout.addWidget(self.target_checkbox, row, 0)
        self.layout.addWidget(self.target_button, row, 1)
        self.layout.addWidget(self.target_filepath, row, 2)
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
        self.layout.addWidget(self.max_missed_edit, row, 1)
        row += 1
        self.layout.addWidget(self.fdr_label, row, 0)
        self.layout.addWidget(self.fdr_edit, row, 1)
        row += 1
        self.layout.addWidget(self.protease_label, row, 0)
        self.layout.addWidget(self.protease_combo_box, row, 1)
        row += 1
        self.layout.addWidget(self.run_button, row, 0)

        self.layout.setAlignment(Qt.AlignTop)
        self.horizontalGroupBox.setAlignment(Qt.AlignTop)

        self.horizontalGroupBox.setLayout(self.layout)


def run_gui():
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())

if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
