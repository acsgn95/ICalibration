"""
THE ICALIBRATION METHOD OF MAGNETOMETERS


The ICalibration method, which was developed jointly by the researchers 
Hairong Yu, Lin Ye, Ying Guo, and Steven Su, was applied.

In this method, the magnetometer is calibrated based on the calibrated 
accelerometer data and the magnetic inclination angle.

Article in IEEE Sensors Journal - May 2020

DOI: 10.1109/JSEN.2020.2995876

Code by : Ahmet Cosgun - Geomatics Engineer

Created : 24 April 2022

"""

import numpy as np
from math import *
import matplotlib.pyplot as plt


class ICalibration():
    
    """
    param :
        
        totalMagneticField : (True total magnetic field)The square root of the sum of the squares of the local magnetic field
        inclinationDegree : Dip angle of magnetik field

        https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm

    
    """

    def __init__(self,totalMagneticField,inclinationDegree):

        self.totalMagneticField = totalMagneticField
        self.inclination = radians(inclinationDegree)
        self.data = None
        self.measMagX = []
        self.measMagY = []
        self.measMagZ = []
        self.accX = []
        self.accY = []
        self.accZ = []
        self.distortion = None
        self.offset = None
        self.X = np.ones((12,9))


    def getData(self,measurements):

        self.data = measurements

        try:
            for meas in self.data:
                self.accX.append(meas[0])
                self.accY.append(meas[1])
                self.accZ.append(meas[2])
                self.measMagX.append(meas[3])
                self.measMagY.append(meas[4])
                self.measMagZ.append(meas[5])

        
        except ValueError:
            print("Review the data! You have to send 6 pieces of data of 12 observations...")

        self.X = np.ones((len(self.data),9))

    def leastSquareMatrix(self):

        for i in range(len(self.data)):

            step = [self.accX[i]*self.measMagX[i],
                    self.accX[i],
                    self.accY[i]*self.measMagY[i],
                    self.accY[i],
                    self.accZ[i]*self.measMagZ[i],
                    self.accZ[i],
                    self.accX[i]*self.measMagY[i]+self.accY[i]*self.measMagX[i],
                    self.accX[i]*self.measMagZ[i]+self.accZ[i]*self.measMagX[i],
                    self.accY[i]*self.measMagZ[i]+self.accZ[i]*self.measMagY[i]]

            for j in range(9):

                self.X[i][j] = step[j]

    
    def calibration(self,data):

        self.getData(data)
        self.leastSquareMatrix()

        L = self.totalMagneticField*cos((pi/2)-self.inclination)
        LMatrix = np.ones((len(self.data),1))*L

        Xtranspose = np.transpose(self.X)

        stepOne = np.dot(Xtranspose,self.X)
        stepTwo = np.linalg.inv(stepOne)
        stepThree = np.dot(stepTwo,Xtranspose)
        beta = np.dot(stepThree,LMatrix)
        self.distortion = np.array([[beta[0][0],beta[6][0],beta[7][0]],
                                    [beta[6][0],beta[2][0],beta[8][0]],
                                    [beta[7][0],beta[8][0],beta[4][0]]])
        
        self.offset = np.array([[beta[1][0]],
                                [beta[3][0]],
                                [beta[5][0]]])
        
        print("\nCALIBRATION PARAMETERS - ICalibration\n")
        print("DISTORTION: \n G: \n",self.distortion,"\n")
        print("OFFSET: \n O: \n",self.offset,"\n")
        print("FOR CALIBRATION: \nmagCalibrated = G*magMeasured + O\n")


