#!/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 10:41:35 2021

@author: aaronsontag

7/28/2023 ACS
    Updating with PIII diagnostic names and path to data folders. Needed to
    update data folder logic to deal with shot numbers less than 6 digits. The
    previous calibration files have been superseded by a new file with more 
    info.
8/9/2023 ACS
    Fixed bugs and added cross-platform machine checks. Cleanup routine now
    mimics P3 igor data display wCleanup routine with matching FL and Bdot
    signals. CalcIp routine needs to be developed for P3 so that routine
    is not updated and just uses plasmarogB.
"""

from os import path
from igor import binarywave
from scipy import signal
from platform import system
import numpy as np

class ShotData:
    def __init__(self, shotnum):
        """
        This loads the data needed for MHD equilibrium reconstruction for a
        specific shot number from the Pegasus Data Archive. The calibration
        file is loaded and the raw data is processed.

        Returns
        -------
        None.

        """
        self.shotnum = shotnum
        self.DataPath = ''

        FluxLoops = ['CFL01', 'CFL02', 'CFL03', 'CFL04', 'CFL05',
                     'CFL06', 'CFL07', 'CFL08', 'CFL09', 'CFL10',
                     'CFL11', 'CFL12', 'CFL13', 'CFL04', 'CFL05',
                     'CTor1', 'CTor2', 'CTor3', 'CTor4', 'CTor5',
                     'CTor6','DiamagA', 'DiamagB', 'DiamagC',
                     'NCFL01', 'NCFL02', 'NCFL03', 'NCFL04',
                     'NCFL05', 'NCFL06', 'NCFL07', 'NCFL08',
                     'NCFL09', 'NCFL10', 'NCFL11', 'NCFL12',
                     'NCFL13', 'NCFL14', 'NCFL15', 'NCFL16',
                     'NCFL17', 'NCFL18', 'NCFL19', 'NCFL20',
                     'WFL01', 'WFL02', 'WFL03', 'WFL04', 'WFL05', 'WFL06',
                     'WFL07', 'WFL08', 'WFL09', 'WFL10', 'WFL11', 'WFL12',
                     'WFL13', 'WFL14', 'WFL15', 'WFL16']

        BDots = ['PDX01', 'PDX02', 'PDX03', 'PDX04', 'PDX05', 'PDX06',
                 'PDX07', 'PDX08', 'PDX09', 'PDX10', 'PDX11', 'PDX12', 
                 'PDX13', 'OTor1', 'OTor2', 'OTor3', 'OTor4', 'OTor5', 
                 'OTor6', 'CTor1', 'CTor2', 'CTor3', 'CTor4', 'CTor5',
                 'CTor6', 'CPA01', 'CPA02', 'CPA03', 'CPA04', 'CPA05', 
                 'CPA06', 'CPA07', 'CPA08', 'CPA09', 'CPA10', 'CPA11', 
                 'CPA12', 'CPA13', 'CPA14', 'CPA15', 'CPA16', 'CPA17', 
                 'CPA18', 'CPA19', 'CPA20', 'CPA21']
        
        Currents = ['PlasmaRogA', 'PlasmaRogB', 'RT_DPWMi_EF123', 
                    'RT_DPWMi_EF45', 'RT_DPWMi_EF678', 'RT_DPWMi_TF']

        self.sigs = {'FluxLoops': FluxLoops, 'BDots': BDots, 'Currents' : Currents}
        self.raw = {'FluxLoops': {}, 'BDots' : {}, 'Currents' : {}}
        self.calData = {'FluxLoops': {}, 'BDots': {}, 'Currents' : {}}
        self.calFile = {}

    def getFolderPath(self):
        """
        shot folder name is expected in PXXXXXX format. This folder is under a
        folder with first two significant digits and a subfolder of the next
        two significant digits.
        For example T109756 is located at DataPathRoot/100000/9700/T109756

        Parameters
        ----------
        shotnum : TYPE
            Shot number.

        Returns
        -------
        None
        """
        # choose data root path based on os type
        machine = system()
        PathDict = {'Darwin' : path.join('/Volumes', 'Pegasus_Data_Archive', 'p3data'),
                    'Windows' : path.join('P:',  'p3data'),
                    'Linux' : path.join('/mnt', 'Pegasus_Data_Archive', 'p3data')}
        DataPathRoot = PathDict[machine]
        
        # this block gives correct folder based on shot number and adds to root path
        shotnum = self.shotnum
        if len(str(shotnum)) == 4:
            shotstr = '00' + str(shotnum)
        elif len(str(shotnum)) == 5:
            shotstr = '0' + str(shotnum)
        else:
            shotstr = str(shotnum)
        f1Str = shotstr[:2] + '0000'
        f2Str = shotstr[2:4] + '00'
        fName = 'T' + shotstr
        PathName = path.join(DataPathRoot, f1Str, f2Str, fName)
        self.DataPath = PathName

    def LoadRaw(self):
        """
        Loads raw data into a data structure. Attempts to load all diagnostic
        names in init class constructor.

        Parameters
        ----------
        shotnum : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        sGroups = self.sigs.keys()
        DataPath = self.DataPath
        print(DataPath)
        for group in sGroups:
            for sig in self.sigs[group]:
                self.raw[group][sig] = IgorLoad(sig, DataPath,quiet=True)

    def LoadCalData(self):
        """
        The signal calibration data is stored in a text file DAS.conf. The 
        beginning of this file has irrelevant header info that needs to be 
        skipped over. Diagnostic calibration information is stored in six
        line blocks that need to be parsed with a blank line in between.

        Returns
        -------
        None.

        """

        calFile = path.join(self.DataPath, 'DAS.conf')

        with open(calFile, 'r', errors='replace') as f:
            DiagLines = f.readlines()[64:]
        DiagNames = DiagLines[0::7]
        CalVals = DiagLines[4::7]
        CalUnits = DiagLines[3::7]
        for i in range(len(DiagNames)):
            self.calFile[DiagNames[i].strip('[]\n')] = [float(CalVals[i].split()[-1]), CalUnits[i].split()[-1]]

    def CalcIp(self):
        """
        This routine takes the raw plasma rogowski signal and calculates Ip.
        For now use the method of integration and subtraction of the
        CoreFluxLoop3 signal. Need a vacuum shot to calibrate background
        subtraction. Will need to be redone for Pegasus III.

        Returns
        -------
        None.

        """
        from numpy import arange, interp

        PlasRogRaw = self.raw['Currents']['PlasmaRogB']
        cleanSig,cleanSigTime = CleanUp(PlasRogRaw, intFlag = True)
        calSig = float(self.calFile['PlasmaRogB'][0]) * cleanSig
        self.calData['Currents']['Ip'] = {}
        self.calData['Currents']['Ip']['data'] = calSig
        self.calData['Currents']['Ip']['time'] = cleanSigTime

    def ProcBdot(self):
        """
        Mirnov processing to get field from raw voltage signal. This requires
        the raw voltage to be multiplied by the value in the calibration file
        as well as by 1/NA from the BdotNA routine. Processed signals are
        stored in the self.calData['Bdots'] tag. Calibrated signals are stored
        with a 'data' key and a 'time' key, which both contain arrays.

        Returns
        -------
        None.

        """

        OneOverNA = BdotNA()
        BDots = self.raw['BDots']
        for sig in BDots.keys():
            rawSig = BDots[sig]
            try:
                cleanSig,cleanSigTime = CleanUp(rawSig, intFlag = True)
                calSig = OneOverNA[sig] * float(self.calFile[sig][0]) * cleanSig
                #pass
            except TypeError:
                print('no data for ' + sig)
                calSig = []
                cleanSigTime = []
            self.calData['BDots'][sig] = {}
            self.calData['BDots'][sig]['data'] = calSig
            self.calData['BDots'][sig]['time'] = cleanSigTime

    def ProcFL(self):
        """
        Flux loop processing. Integrating and signal cleanup. This requires
        the raw voltage to be multiplied by the value in the calibration file.
        Processed signals are stored in the self.calData['Bdots'] tag. 
        Calibrated signals are stored with a 'data' key and a 'time' key, 
        which both contain arrays.

        Returns
        -------
        None.

        """

        FluxLoops = self.raw['FluxLoops']
        for sig in FluxLoops.keys():
            rawSig = FluxLoops[sig]
            try:
                cleanSig,cleanSigTime = CleanUp(rawSig, intFlag = True)
                calSig = float(self.calFile[sig][0]) * cleanSig
            except TypeError:
                print('no data for ' + sig)
                calSig = []
                cleanSigTime = []
            self.calData['FluxLoops'][sig] = {}
            self.calData['FluxLoops'][sig]['data'] = calSig
            self.calData['FluxLoops'][sig]['time'] = cleanSigTime
            
    def TestPlot(self):
        """
        Simple plotting routine to spot check data.

        Returns
        -------
        None.

        """
        
        plotx = self.calData['Currents']['Ip']['time']
        ploty = self.calData['Currents']['Ip']['data']
            
        from matplotlib.pyplot import plot, show
        plot(plotx, ploty)
        show()
        
    def LoadProcData(self):
        """
        Exports processed data structure

        Returns
        -------
        None.

        """
        return self.calData


def IgorLoad(signame, pathname, quiet=True):
    """
    Given a signame which is an IgorPro .ibw save file, this routine will
    load the stored data and return a dictionary containing the deltaX and
    Y values.

    Parameters
    ----------
    signame : TYPE
        DESCRIPTION.
    pathname : TYPE
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    sig = signame + '.ibw'
    try:
        data = binarywave.load(path.join(pathname, sig))
        sname = data['wave']['wave_header']['bname']
        try:
            deltaX = data['wave']['wave_header']['hsA']
        except KeyError:
            deltaX = data['wave']['wave_header']['sfA'][0]
            #print(sname, deltaX)
        wData = data['wave']['wData']
        if not quiet:
            print(sig + ' loaded  deltaX = ' + str(deltaX))
        return {'deltaX': deltaX, 'data': wData}
    except FileNotFoundError:
        print(sig + ' not loaded - file not found')
        return {'deltaX': [], 'data': []}
    

def CleanUp(signal, intFlag = False):
    """
    This procedure performs the same tasks as the Pegasus Igor Data Loader 
    wCleanup function. The input signal is in the form loaded from igor pro
    with an array of y values and a delta x. The processed signal will be 
    output as an array of (x,y) pairs.

    Parameters
    ----------
    signal : raw signal loaded from igor containg a data key that has the y
        values and a deltaX key that has the delta x for the data spacing.
    
    intFlag : boolian flag to digitally integrate signal

    Returns
    -------
    CleanData : (x,y) array of processed data
    An array of the data with offset and drift removed. Optional integration
    if desired.

    """
    from scipy import integrate
    
    # turn igor y & deltaX representation into x & y arrays
    sig = signal['data']
    dx = signal['deltaX']
    npts = len(sig)
    sigTime = dx * np.arange(npts)

    # subtract average over first 1 ms to remove baseline offset
    dStart = 0.0
    dEnd = 0.001
    dRange = np.flatnonzero((sigTime > dStart) & (sigTime < dEnd))
    offSig = sig - np.mean(sig[dRange])
    
    # integrate if desired
    if intFlag:
        intSig = integrate.cumtrapz(offSig, x=sigTime, initial=0)
        offSig = intSig
    
    # remove linear baseline drift - subtracting line from whole array
    cleanSig = offSig
    for i in range(0,len(offSig)): 
        cleanSig[i] -= (offSig[-1] - offSig[0]) * sigTime[i] / (sigTime[-1] - sigTime[0])
    
    return cleanSig,sigTime

def BdotNA():
    """
    The BdotDict gives the values for 1/NA based on probe names. The 
    values have been updated to reflect the P3 calibrations as of 
    7/31/2023. - ACS

    Parameters
    ----------

    Returns
    -------
    BdotDict : Dictionary with Pegasus III B-dot coil NA calibration data

    """
    BdotDict = {"PDX01": 31.40986856, "PDX02": 30.61780031,
               "PDX03": 30.83447292, "PDX04": 30.18429813,
               "PDX05": 29.80654102, "PDX06": 31.5045286,
               "PDX07": 29.94459065, "PDX08": 30.03320686,
               "PDX09": 30.94971522, "PDX10": 31.56406303,
               "PDX11": 30.99254654, "PDX12": 29.48266806,
               "PDX13": 31.41396525, "OTor1": 51.43020212,
               "OTor2": 53.77402821, "OTor3": 53.68001112,
               "OTor4": 50.42248514, "OTor5": 52.05609586,
               "OTor6": 51.94900848, "CTor1": 176.3906218,
               "CTor2": 202.8254814, "CTor3": 183.3700253,
               "CTor4": 185.2652264, "CTor5": 187.4727088,
               "CTor6": 214.1072423, 
               "CPA01": 179.71749, "CPA02": 198.80251,
               "CPA03": 197.43882, "CPA04": 185.42609,
               "CPA05": 186.73373, "CPA06": 190.3981965,
               "CPA07": 192.9140393, "CPA08": 198.2726724,
               "CPA09": 182.9567595, "CPA10": 193.7743813,
               "CPA11": 196.2893948, "CPA12": 185.3693143,
               "CPA13": 191.0163613, "CPA14": 192.7280717,
               "CPA15": 172.6144269, "CPA16": 171.9967657,
               "CPA17": 182.0353776, "CPA18": 209.1828164,
               "CPA19": 182.0272987, "CPA20": 184.948528,
               "CPA21": 174.5560926}

    return BdotDict

def main(shot):
    """
    Main program for command line calls

    Parameters
    ----------
    shot : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    current = ShotData(shot)
    current.getFolderPath()
    current.LoadRaw()
    current.LoadCalData()
    current.CalcIp()
    current.ProcFL()
    current.ProcBdot()
    current.TestPlot()

# --- Launch main() ---------------------------------------------------------
if __name__ == '__main__':
    import argparse
    #import configargparse.ArgumentParser as AP
    #import configargparse.RawDescriptionHelpFormatter as HF

    parser = argparse.ArgumentParser(description='loads/processes Pegasus data',
                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('shot', help="shotnumber", type=int)

    args = parser.parse_args()

    shot = args.shot

    main(shot)
