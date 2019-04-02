#####################

#from fbs_runtime.application_context import ApplicationContext
#from PyQt5.QtWidgets import QMainWindow

from __future__ import division
from PyQt5 import QtCore, QtWidgets, QtGui

from biota_gui import Ui_MyPluginDialogBase  # importing our generated file
import sys
import os
import platform



#####
# A bunch of functions

def install_biota(directory):
    if os.path.isdir(directory + 'biota.egg-info/') == False:
        install_msg()
        os.system("cd " + directory + " && python setup.py install")


def install_msg():
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Information)
    msg.setText("I didn't find Biota")
    msg.setInformativeText("Let me install it for you")
    msg.setWindowTitle("Error")
    msg.exec_()


def download_data(directory, Lat, Lon, y):
    # Make the file name
    if int(Lat)>0: Lat_desc = 'N'
    else: Lat_desc = 'S'
    if int(Lon)>0: Lon_desc = 'E'
    else: Lon_desc = 'W'
    Lat_str = str(abs(int(Lat))).zfill(2); Lon_str = str(abs(int(Lon))).zfill(3); Y_str = y[-2:]
    folder = Lat_desc + Lat_str + Lon_desc + Lon_str + '_' + Y_str

    download = True
    for root, dirs, files in os.walk(directory):
        if root[:-4] == directory + folder:
            download = False

    return download




def error_msg(name_str, info_str):
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.setText(name_str)
    msg.setInformativeText(info_str)
    msg.setWindowTitle("Error")
    msg.exec_()


def success_msg():
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Information)
    msg.setText("I'm all done!")
    msg.setInformativeText("Have a nice day :-) ")
    msg.setWindowTitle("Success")
    msg.exec_()


def is_digit(n):
    try:
        int(n)
        return True
    except ValueError:
        return False



def save_params(INP, OUT,LAT, LON, SIZE, bool_F_prop, bool_F_chge, Y1, Y2, bool_Gammao, bool_AGB, bool_WCV, bool_BCH, bool_FCH, bool_RSK, bool_FIL, SUB, POL, ATH, BTH, ACH, BCH, PCH):
    file = open(OUT+"parameters.txt","w")

    file.write("Input_folder = " + INP + " \n")
    file.write("Output_folder = " + OUT + " \n")
    file.write("Latitude = " + LAT + " deg \n")
    file.write("Longitude = " + LON + " deg \n")
    file.write("Tile_size = " + SIZE + " deg \n")
    file.write("Forest_properties = " + str(bool_F_prop) + " \n")
    file.write("Forest_change = " + str(bool_F_chge) + " \n")
    file.write("Year 1 = " + Y1 + " \n")
    file.write("Year 2 = " + Y2 + " \n")
    file.write("Gamma0 = " + str(bool_Gammao) + " \n")
    file.write("AGB = " + str(bool_AGB) + " \n")
    file.write("Forest_cover = " + str(bool_WCV) + " \n")
    file.write("Biomass_change = " + str(bool_BCH) + " \n")
    file.write("Forest_cover_change = " + str(bool_FCH) + " \n")
    file.write("Deforestation_risk = " + str(bool_RSK) + " \n")
    file.write("Filter = " + str(bool_FIL) + " \n")
    file.write("Resample_factor = " + SUB + " \n")
    file.write("Polarisation = " + POL + " \n")
    file.write("Area_threshold = " + ATH + " \n")
    file.write("Biomass_threshold = " + BTH + " \n")
    file.write("Area_change_threshold = " + ACH + " \n")
    file.write("Biomass_change_treshold = " + BCH + " \n")
    file.write("Percentage_change_threshold = " + PCH + " \n")

    file.close()



#####

class mywindow(QtWidgets.QMainWindow):

    def __init__(self):
        super(mywindow, self).__init__()
        self.ui = Ui_MyPluginDialogBase()
        self.ui.setupUi(self)

        # Set up the text wrap
        self.ui.label_22.setWordWrap(True)
        self.ui.label_23.setWordWrap(True)

        # Set up the input folder button
        self.ui.lineEdit_INP.clear()
        self.ui.pushButton_INP.clicked.connect(self.get_input_dir)

        # Set up the ouput folder button
        self.ui.lineEdit_OUT.clear()
        self.ui.pushButton_OUT.clicked.connect(self.get_output_dir)

        # Set up the Map button
        self.ui.pushButton_3.clicked.connect(self.maplink)

        # Set up the Map button
        self.ui.helpButton.clicked.connect(self.helplink)


        # Set-up default data location
        self.ui.lineEdit_INP.setText(os.getcwd() + '/DATA/')
        self.ui.lineEdit_OUT.setText(os.getcwd() + '/DATA/')


        # Set up the logos
        SMFM_label = self.ui.label_21
        SMFM_pixmap = QtGui.QPixmap('gui/Logos/SMFM_Logo.png')
        SMFM_label.setScaledContents(True); SMFM_label.setPixmap(SMFM_pixmap)

        UoE_label = self.ui.label_28
        UoE_pixmap = QtGui.QPixmap('gui/Logos/UoE_Logo.jpg')
        UoE_label.setScaledContents(True); UoE_label.setPixmap(UoE_pixmap)

        LTS_label = self.ui.label_24
        LTS_pixmap = QtGui.QPixmap('gui/Logos/LTS_Logo.png')
        LTS_label.setScaledContents(True); LTS_label.setPixmap(LTS_pixmap)

        white_label = self.ui.label_25
        white_pixmap = QtGui.QPixmap('gui/Logos/white.jpg')
        white_label.setScaledContents(True); white_label.setPixmap(white_pixmap)

        white_label = self.ui.label_26
        white_pixmap = QtGui.QPixmap('gui/Logos/white.jpg')
        white_label.setScaledContents(True); white_label.setPixmap(white_pixmap)


        # Set-up the help buttons
        self.ui.helpButton_lat.clicked.connect(self.lat_help)
        self.ui.helpButton_lon.clicked.connect(self.lon_help)
        self.ui.helpButton_g0.clicked.connect(self.g0_help)
        self.ui.helpButton_bio.clicked.connect(self.bio_help)
        self.ui.helpButton_for.clicked.connect(self.for_help)
        self.ui.helpButton_bch.clicked.connect(self.bch_help)
        self.ui.helpButton_fch.clicked.connect(self.fch_help)
        self.ui.helpButton_rsk.clicked.connect(self.rsk_help)
        self.ui.helpButton_fil.clicked.connect(self.fil_help)
        self.ui.helpButton_sub.clicked.connect(self.sub_help)
        self.ui.helpButton_pol.clicked.connect(self.pol_help)
        self.ui.helpButton_ath.clicked.connect(self.ath_help)
        self.ui.helpButton_bth.clicked.connect(self.bth_help)
        self.ui.helpButton_ach.clicked.connect(self.ach_help)
        self.ui.helpButton_BCH.clicked.connect(self.BCH_help)
        self.ui.helpButton_pch.clicked.connect(self.pch_help)





    def get_input_dir(self):
        filename = QtWidgets.QFileDialog.getExistingDirectory(self, "Select input folder ")
        filename = filename + "/"
        self.ui.lineEdit_INP.setText(filename)
        self.ui.lineEdit_OUT.setText(filename)



    def get_output_dir(self):
        filename = QtWidgets.QFileDialog.getExistingDirectory(self, "Select input folder ")
        filename = filename + "/"
        self.ui.lineEdit_OUT.setText(filename)



    def maplink(self):
        type = 'data=!3m1!1e3'
        lat = self.ui.lineEdit_Lat.text()
        lon = self.ui.lineEdit_Lon.text()
        if self.ui.horizontalSlider.value() == 1: scale = '8292m'
        else: scale = '46833m'

        if [lat,lon] != ['','']:
            QtGui.QDesktopServices.openUrl(QtCore.QUrl('https://www.google.co.uk/maps/@'+lat+','+lon+','+scale+'/'+type))
        else:
            QtGui.QDesktopServices.openUrl(QtCore.QUrl('https://www.google.co.uk/maps/'))



    def helplink(self):
        QtGui.QDesktopServices.openUrl(QtCore.QUrl('https://biota.readthedocs.io/en/latest/'))


    # help functions
    def lat_help(self):
        self.ui.label_22.setText('Angular distance North (+) or South (-) of the Equator. Please enter the latitude of the top left corner of your area of interest.')
    def lon_help(self):
        self.ui.label_22.setText('Angular distance East (+) or West (-) of the Greenwich Meridian. Please enter the longitude of the top left corner of your area of interest.')
    def g0_help(self):
        self.ui.label_22.setText('Radar backscatter from ALOS data. Biomass is calculated from horizontal-send and horizontal-receive (HV) polarisation')
    def bio_help(self):
        self.ui.label_22.setText('Above-ground biomass calculated by BIOTA using HV polarisation. The calibration is ??????????')
    def for_help(self):
        self.ui.label_22.setText('Forest cover (binary result) determined by combining the above-ground biomass with the area threshold and biomass threshold for forested areas.')
    def bch_help(self):
        self.ui.label_22.setText('Change in above-ground biomass between year 1 and year 2')
    def fch_help(self):
        self.ui.label_22.setText('Change in forest cover between year 1 and year 2, determined by combining change in biomass with the area/biomass/percentage change parameters. Area threshold and biomass threshold also influence the result by affecting the calcuation of forest cover for each year.')
    def rsk_help(self):
        self.ui.label_22.setText('Probability of deforestation based on ?????????????')
    def fil_help(self):
        self.ui.label_22.setText('Speckle filter to get rid of confusing white specks on your final image.')
    def sub_help(self):
        self.ui.label_22.setText('Divides the final raster resolution by the integer in the box.')
    def pol_help(self):
        self.ui.label_22.setText('Four polarisation options. Polarisation is the process through which rays are sent and received with a certain orientation (horizontal or vertical in this case). HV polarisation is used to calculate above-ground biomass.')
    def ath_help(self):
        self.ui.label_22.setText('Minimum area for the forest cover to consider a patch of pixels to be forested.')
    def bth_help(self):
        self.ui.label_22.setText('Minimum above-ground biomass density for the forest cover to consider a patch of pixels to be forested.')
    def ach_help(self):
        self.ui.label_22.setText('Minimum change in forested area for the change algorithm to classify biomass loss as deforestation.')
    def BCH_help(self):
        self.ui.label_22.setText('Minimum change in forested biomass for the change algorithm to classify biomass loss as deforestation.')
    def pch_help(self):
        self.ui.label_22.setText('Minimum relative change in forested biomass for the change algorithm to classify biomass loss as deforestation.')



    def reject(self):
        self.close()


    def accept(self):
        # Set up the paths
        biota_inst = os.getcwd() + '/'
        biopath = os.getcwd() + '/cli/'

        # Install biota if necessary
        install_biota(biota_inst)

        # Reset the bars
        self.ui.progressBar_DLD.setValue(0)
        self.ui.progressBar_PRC.setValue(0)


        # Get the parameters from the lineEdits and the spinBox
        INP = str(self.ui.lineEdit_INP.text())
        OUT = str(self.ui.lineEdit_OUT.text())

        LAT = str(self.ui.lineEdit_Lat.text())
        LON = str(self.ui.lineEdit_Lon.text())
        SIZE = str(self.ui.horizontalSlider.value())
        Y1 = str(self.ui.lineEdit_Y1.text())
        Y2 = str(self.ui.lineEdit_Y2.text())

        SUB = str(self.ui.spinBox_SUB.text())

        POL = str(self.ui.lineEdit_POL.text())
        PCH = str(self.ui.lineEdit_PCH.text())
        BTH = str(self.ui.lineEdit_BTH.text())
        BCH = str(self.ui.lineEdit_BCH.text())
        ATH = str(self.ui.lineEdit_ATH.text())
        ACH = str(self.ui.lineEdit_ACH.text())


        # Get the bools from the checkBoxes
        bool_F_prop = self.ui.checkBox_pty.isChecked()
        bool_F_chge = self.ui.checkBox_chg.isChecked()
        bool_Gammao = self.ui.checkBox_Go.isChecked()
        bool_AGB = self.ui.checkBox_AGB.isChecked()
        bool_FCH = self.ui.checkBox_FCH.isChecked()
        bool_WCV = self.ui.checkBox_WCV.isChecked()
        bool_BCH = self.ui.checkBox_BCH.isChecked()
        bool_RSK = self.ui.checkBox_RSK.isChecked()

        bool_FIL = self.ui.checkBox_FIL.isChecked()



        # Check the presence/format of all required inputs and print error messages

        good2run = True

        if INP == "" or os.path.isdir(INP) == False :
            error_msg('Input Directory Error', 'Please enter a valid input directory'); good2run = False
        elif OUT == "" or os.path.isdir(INP) == False:
            error_msg('Output Directory Error', 'Please enter a valid output directory'); good2run = False
        elif is_digit(LAT) == False or is_digit(LON) == False:
            error_msg('Lat/Lon Error', 'Please tell me where to work'); good2run = False
        elif bool_F_prop == False and bool_F_chge == False:
            error_msg('Output Error', 'I have nothing to do'); good2run = False
        elif bool_F_prop == True:
            if Y1.isdigit() == False :
                error_msg('Year 1 Error', 'Please enter a valid year'); good2run = False
            if [bool_Gammao, bool_AGB, bool_WCV] == [False, False, False]:
                error_msg('Output Error', 'Please choose an output'); good2run = False




        if bool_F_chge == True:
            if Y2.isdigit() == False or int(Y2) <= int(Y1):
                error_msg('Year 2 Error', 'Please choose a valid second year'); good2run = False
            elif [bool_BCH, bool_FCH, bool_RSK] == [False, False, False]:
                error_msg('Output Error', 'Please choose an output'); good2run = False




        # Check the validity of optional inputs and restore default if necessary
        if any(POL in s for s in ['HV', 'VV', 'VH', 'HH']) == False:
            self.ui.lineEdit_POL.setText('HV'); POL = str(self.ui.lineEdit_POL.text())
        if BTH.isdigit() == False:
            self.ui.lineEdit_BTH.setText('10'); BTH = str(self.ui.lineEdit_BTH.text())
        if BCH.isdigit() == False:
            self.ui.lineEdit_BCH.setText('0'); BCH = str(self.ui.lineEdit_BCH.text())
        if ACH.isdigit() == False:
            self.ui.lineEdit_ACH.setText('0'); ACH = str(self.ui.lineEdit_ACH.text())
        if ATH.isdigit() == False:
            self.ui.lineEdit_ATH.setText('1'); ATH = str(self.ui.lineEdit_ATH.text())
        if PCH.isdigit() == False:
            self.ui.lineEdit_PCH.setText('0'); PCH = str(self.ui.lineEdit_PCH.text())


        # Get the filtering parameter
        if bool_FIL == True: FIL = ""
        else: FIL = " -nf "

        # Get the tile size
        if SIZE == "1": L = ""
        else: L = " -l "


        # Get the years list
        Y_list = [Y1]
        if Y2.isdigit(): Y_list.append(Y2)

        # Get the output list
        Out_list = []
        if bool_Gammao == True: Out_list.append('Gamma0')
        if bool_AGB == True: Out_list.append('AGB')
        if bool_WCV == True: Out_list.append('WoodyCover')

        # Get the change list
        Chg_list = []
        if bool_BCH == True: Chg_list.append('AGBChange')
        if bool_FCH == True: Chg_list.append('ChangeType')
        if bool_RSK == True: Chg_list.append('RiskMap')


        # Run the code
        # (for now, make sure you have the biota environment activated)
        print ('Are we good to go?', good2run)

        if good2run == True:

            # download the data if needed
            dld_count = 0.
            for y in Y_list:
                if download_data(INP, LAT, LON, y) == True:
                    os.system("python " + biopath + "download.py -lat " + LAT + " -lon " + LON + L + " -y " + y + " -o " + INP + " -r")
                    dld_count +=1.
                    self.ui.progressBar_DLD.setValue(int(100*dld_count/(1.0*len(Y_list))))
            self.ui.progressBar_DLD.setValue(100)

            # run the properties
            prc_count = 0.
            if bool_F_prop == True:
                self.ui.label_12.setText('Processing data (property)')
                for y in Y_list:
                    for o in Out_list:
                        os.system("python " + biopath + "property.py -dir " + INP + " -lat " + LAT + " -lon " + LON + " -y " + y + " -o " + o + " " + FIL + " -ds " + SUB + " -od " + OUT + " -p " + POL + " -ft " + BTH + " -at " + ATH )
                        prc_count +=1.
                        self.ui.progressBar_PRC.setValue(int(100*prc_count/ (1.0*len(Y_list)*len(Out_list))))

            # run the change
            prc_count = 0.
            if bool_F_chge == True:
                self.ui.label_12.setText('Processing data (change)')
                self.ui.progressBar_PRC.setValue(0)
                for o in Chg_list:
                    os.system("python " + biopath + "change.py -dir " + INP + " -lat " + LAT + " -lon " + LON + " -y1 " + Y1 + " -y2 " + Y2 + " -o " + o + " " + FIL + " -ds " + SUB + " -od " + OUT + " -ft " + BTH + " -at " + ATH + " -ct " + ACH + " -mt " + BCH + " -it " + PCH )
                    prc_count +=1.
                    self.ui.progressBar_PRC.setValue(int(100*prc_count/ (1.0*len(Chg_list))))


            self.ui.progressBar_PRC.setValue(100)


            # Print the input parameters in a paramfile
            save_params(INP, OUT,LAT, LON, SIZE, bool_F_prop, bool_F_chge, Y1, Y2, bool_Gammao, bool_AGB, bool_WCV, bool_BCH, bool_FCH, bool_RSK, bool_FIL, SUB, POL, ATH, BTH, ACH, BCH, PCH)

            success_msg()









# KEEP THIS BIT!

if __name__ == "__main__":

    app = QtWidgets.QApplication([])
    win = mywindow()
    layout = QtWidgets.QVBoxLayout()

    app.setStyle('Fusion')

    win.show()
    app.exec_()
