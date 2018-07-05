# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 09:46:22 2018

@author: hatlab
"""

#for 0IF tuning, 


import numpy as np
from pulse import Gaussian, Square, Marker
from sequence import Sequence
from timeit import default_timer
import AWGFile
import qt
import AWGclean as clean


def setAWGpalse( set_ch3offset = 0.02, set_ch4_offset = -0.002, set_ssb_freq = 0.0, set_iqscale = 0.68/1.0, set_phase = np.pi/2, set_skewphase = (0.25)*2*np.pi ): # Ch4/Ch3
    
    time1 = default_timer()
    
    
    ### Parameters of Square wave
    square_widthQ = 1000
    square_heightQ = 5000       # non-zeor value for qubit drive, then with the non-zero ssb, one can tune the qubit drive mixer
    
    square_widthC = 1000
    square_heightC = 0          # 0 DAC value for cavity drive, so that one can tune the DC leakage for the cavity
    
    
    ### Parameters of Markers
    # Marker Q controls the qubit swtich which lets qubits pulses through when on 
    markerQ_delay = 10
    markerQ_width = square_widthQ + 2*markerQ_delay
    markerQ_on = 0#5
    markerQ_off = markerQ_width 
    # Marker C controls the cavity swtich which lets cavity pulses through when on 
    markerC_delay = 10 
    markerC_width = square_widthC + 2*markerC_delay
    markerC_on = 0
    markerC_off = markerC_width #- 5
    #This marker triggers the Alazar to begin aquiring data 
    alazar_marker_width = 100
    alazar_marker_on = 5
    alazar_marker_off = 95
    ###
    
    ### Pulse init
    
    # Parameters
    ssb_freq = 0.0
    iqscale = 0.68/1.0 # Ch4/Ch3
    phase = np.pi/2
    skewphase = (0.25)*2*np.pi
    #
    
    ###############################################################################
    
    # This following block is to create all the pulses that we are going to use 
    # in the measurement. 
    # Note that two pulse with same pulse type but different pulse parameters should 
    # be considered as differet pulses. (e.g Two Gaussian with different amp is two pulses)
    
    pulsenum = 6
    pulse = np.empty(pulsenum, dtype = object)
    i = 0
    while i < pulsenum:
        # This is pulse[1], the square pulse for qubit
        pulse[i] = Square(i, square_widthQ, ssb_freq, iqscale, phase, square_heightQ, skew_phase = skewphase)
        pulse[i].data_generator()
        i = i + 1
        # This is pulse[2], the square pulse for the cavity 
        pulse[i] = Square(i, square_widthC, 0.0, 1, 0, square_heightC)  # cavity drive does not need ssb
        pulse[i].data_generator()    
        i = i + 1
        # This is the marker for the qubit
        pulse[i] = Marker(i, markerQ_width, 5, markerQ_on, markerQ_off)
        pulse[i].data_generator() 
        i = i + 1
        # This is the marker for the cavity
        pulse[i] = Marker(i, markerC_width, 5, markerC_on, markerC_off)
        pulse[i].data_generator() 
        i = i + 1
        # This is the trigger for the alazar 
        pulse[i] = Marker(i, alazar_marker_width, 6, alazar_marker_on, alazar_marker_off)
        pulse[i].data_generator() 
        i = i + 1
        pulse[i] = Marker(i, markerQ_width, 2, markerQ_on, markerQ_off)
        pulse[i].data_generator() 
        i = i + 1
    ###############################################################################
    
    
    ###############################################################################
    ### Parameters for the sequence
    start = 0
    wait_time = 500
    ###
    
    
    ### Sequence 
    shotsnum = 1 # This is how many shots we will do in this experiment, i.e how many wait triggers you need in one run 
    shotsnum = shotsnum + 1 # TODO: there is something wrong with the code, the shots number always have one shift...
    totwavenum = 7*shotsnum # This is how many waves are going to be used in the whole sequence
    sequence = Sequence(shotsnum, totwavenum, pulse)
    wait_offset = 2
    # The following while loop is where you create you sequence. 
    # Tell the sequence where you want a certain pulse to start at which channel
    
    i = 0
    while i < shotsnum:
        sequence.get_block(pulse[5].name, start, channel = 1, wait_trigger = True)
        sequence.get_block(pulse[2].name, start, channel = 1)    
    #    sequence.get_block(pulse[4].name, start, channel = 3)
        
        
        start = start + markerQ_delay
        sequence.get_block(pulse[0].name, start, channel = 3)
        sequence.get_block(pulse[1].name, start, channel = 3)
        sequence.get_block(pulse[2].name, start - markerQ_delay, channel = 3)
        start = start + pulse[0].width
       
        i = i + 1
    ###
    
    ### Call difference method from Sequence to do the calculation to find where to put the waveforms
    i = 0 
    while i < shotsnum-1:
        sequence.sequence_block[i].make_block()
        i = i + 1
    
    AWG = AWGFile.AWGFile(filename)
    sequence.sequence_upload(AWG)
    ####
    
    
    AWGInst.restore(awgname)
    print 'Confucius says it takes the following seconds to finish the code'
    print default_timer() - time1
    
    AWGInst.channel_on(1)
    AWGInst.channel_on(2)
    AWGInst.channel_on(3)
    AWGInst.channel_on(4)
    change_setting = True
    
    if change_setting:
        # Set the DC offset
        ch1_offset = -0.084
        ch2_offset = -0.175
        ch3_offset = set_ch3_offset
        ch4_offset = set_ch4_offset
        
        AWGInst.set_ch1offset(ch1_offset)
        AWGInst.set_ch2offset(ch2_offset)
        AWGInst.set_ch3offset(ch3_offset)
        AWGInst.set_ch4offset(ch4_offset)
        
        # Set the channel voltage
        ch1_amp = 1.0
        ch2_amp = 1.0
        ch3_amp = 1.0
        ch4_amp = 1.0
        
        AWGInst.set_ch1amp(ch1_amp)
        AWGInst.set_ch2amp(ch2_amp)
        AWGInst.set_ch3amp(ch3_amp)
        AWGInst.set_ch4amp(ch4_amp)
    
    AWGInst.run()
    

#==============================================================================
# channel_offset(self, channel_num, offset)
# channel_skew(self, channel_num, skew)
# channel_amp(self, channel_num, amp)
#==============================================================================

    
def findmin( channel_num, sweep_type = 'offset', wait_time = 0.05, sweepoints ):
    #general findmin for all channel & type & sweeppoints
    #in this step we have already set the marker frequency
    peak_value = np.zeros(len(sweeepoints))
    for i in range(len(sweepoints)):
        if sweep_type == 'offset':
            AWGInst.channel_offset(channel_num, sweepoints[i])
        elif sweep_type == 'amp'
            AWGInst.channel_amp(channel_num, sweepoints[i])
        elif sweep_type == 'skew'
            AWGInst.channel_skew(channel_num, sweepoints[i])
            
        qt.msleep(wait_time)
        peak_value[i] = MXA.marker_Y_value(1)

    return sweepoints[np.argmin(peak_value)]

def findclosest( channel_num, sweep_type = 'offset', wait_time = 0.05, sweepoints, ref ):
    #general findclosest for all channel & type & sweeppoints
    #in this step we have already set the marker frequency
    peak_value = np.zeros(len(sweeepoints))
    for i in range(len(sweepoints)):
        if sweep_type == 'offset':
            AWGInst.channel_offset(channel_num, sweepoints[i])
        elif sweep_type == 'amp'
            AWGInst.channel_amp(channel_num, sweepoints[i])
        elif sweep_type == 'skew'
            AWGInst.channel_skew(channel_num, sweepoints[i])
            
        qt.msleep(wait_time)
        peak_value[i] = abs( MXA.marker_Y_value(1) - ref )

    return sweepoints[np.argmin(peak_value)]


#==============================================================================
# tune up DC-offset
# while...
# channel4
# channel3    
# record offset and set AWG sequence.   
#==============================================================================
   
def DC_offset_tune_up( I_channel_num, I_offset, Q_channel_num, Q_offset ):
    #need to set start offset    
    for avgtimes, step in [ (1, 1), (2, 0.1), (5, 0.02), (20, 0.005) ]:  #, (100, 0.001)       
        MXA.set_max_count(avgtimes)
        qt.msleep( avgtimes * 0.05 )
        MXA.new_peak(1)
        I_offset = findmin( I_channel_num, avgtimes * 0.05, np.liespace( I_offset - step, I_offset + step, 10 ) )
        Q_offset = findmin( Q_channel_num, avgtimes * 0.05, np.liespace( Q_offset - step, Q_offset + step, 10 ) )
    
    return I_offset, Q_offset
    
#==============================================================================
# tune up iqscale
# set phase to 0 or pi/2
# channel3
# record iqscale and set sequence
#==============================================================================
   
   
def IQscale_tune_up( I_channel_num, Q_channel_num, IQ_scale, ref ):
    #need to set start IQscale & referance( only I peak level)
    for avgtimes, step in [ (1, 1), (2, 0.1), (5, 0.02), (20, 0.005) ]:         
        MXA.set_max_count(avgtimes)
        qt.msleep( avgtimes * 0.05 )
        MXA.new_peak(1)
        IQ_scale = findclosest( Q_channel_num, 'amp', avgtimes * 0.05, np.liespace( IQscale - step, IQscale + step, 10 ), ref )
    
    return IQscale
   
   
#==============================================================================
# tuneup skewphase
# set phase to pi/4
# channel3
# record skewphase and set sequence
#==============================================================================
   
def Skew_tune_up( I_channel_num, Q_channel_num, Skew, ref ):
    #need to set start IQscale & referance( only I peak level)
    for avgtimes, step in [ (1, 1), (2, 0.1), (5, 0.02), (20, 0.005) ]:         
        MXA.set_max_count(avgtimes)
        qt.msleep( avgtimes * 0.05 )
        MXA.new_peak(1)
        Skew = findclosest( Q_channel_num, 'amp', avgtimes * 0.05, np.liespace( ( Skew - step ) * 2 * np.pi, ( Skew + step ) * 2 * np.pi, 10 ), ref )
    
    return Skew
    
    
    
#==============================================================================
# repeat step 1 & 2 & 3
#==============================================================================

### Define the file name and instrument
name = 'Leakage_check.Awg'
filename = 'C:\\Users\\Public\\' + name
awgname = '\\\HATLAB_3-PC\\Users\\Public\\' + name
AWGInst = qt.instruments.get('AWG')

MXA = qt.instruments.get('MXA')
MXA.set_frequency_center(7.4e9)
MXA.set_frequency_span(20e6)
###

clean.clean()

Ichannel = 4, Qchannel = 3
Ioffset = 0.02, Qoffset = -0.002
IQscale = 1.0
Skew = 0.0
Iterations = 3

for i in range(Iterations):
    setAWGpalse( set_ch3offset = Qoffset, set_ch4_offset = Ioffset, set_iqscale = IQscale, set_phase = np.pi/2, set_skewphase = Skew * 2 * np.pi )
    Ioffset, Qoffset = DC_offset_tune_up( I_channel_num, I_offset, Q_channel_num, Q_offset ):
    
    setAWGpalse(  set_ch3offset = Qoffset, set_ch4_offset = Ioffset, set_iqscale = IQscale, set_phase = np.pi/2, set_skewphase = Skew * 2 * np.pi )
    MXA.set_max_count(100)
    qt.msleep(5)
    MXA.new_peak(1)
    ref = MXA.marker_Y_value(1)
    setAWGpalse( set_ch3offset = Qoffset, set_ch4_offset = Ioffset, set_iqscale = IQscale, set_phase = 0, set_skewphase = Skew * 2 * np.pi )
    IQscale = IQscale_tune_up( I_channel_num, Q_channel_num, IQscale, ref )
    
    setAWGpalse( set_ch3offset = Qoffset, set_ch4_offset = Ioffset, set_iqscale = IQscale, set_phase = np.pi/4, set_skewphase = Skew * 2 * np.pi )
    Skew = Skew_tune_up( I_channel_num, Q_channel_num, Skew, ref )

print Ichannel, Qchannel
print Ioffset, Qoffset
print IQscale
print Skew

   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    
