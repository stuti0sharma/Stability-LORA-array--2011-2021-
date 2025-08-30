import numpy as np
from datetime import datetime, timedelta
from datetime import date
from optparse import OptionParser
import os
import sys
sys.path.append('//vol/astro7/lofar/kmulrey/software/root/root_install/lib') 
import ROOT
import LORAparameters as LORA
import os.path
from datetime import date
import detector
import event


nTraceV1=4000
nTraceV2=0
nDetV1=20
dtV1=2.5
timesV1=np.arange(0,4000*dtV1,dtV1)

avg_thresholdV1=np.asarray([ 27.24867982,  27.08347496,  21.3469044,   22.45092564,  16.96505311,
16.62127466,  17.87311077,  17.90886191,  20.3131563,   20.87103187,
20.18339909,  22.10192716,  22.96687405,  23.04840668,  13.57128983,
13.92050076,  27.7231563,   15.69400607,  23.10057663,  24.20940819, 
20.76010091, 20.76010091, 20.76010091, 20.76010091, 20.76010091, 
20.76010091, 20.76010091, 20.76010091, 20.76010091, 20.76010091,
20.76010091, 20.76010091, 20.76010091, 20.76010091, 20.76010091,
20.76010091, 20.76010091, 20.76010091, 20.76010091, 20.76010091]) # The last 20 threshold were missing, they are now the mean of the other thresholds combined. Awaiting a measured avg threshold


def all_events_in_root(tag, data_dir):

    timestamps = []
    ns_timestamps = []

    timestamps = np.loadtxt(data_dir+tag+'_array.log',usecols = (0))
    ns_timestamps = np.loadtxt(data_dir+tag+'_array.log',usecols = (1))
    n_actived_detectors = np.loadtxt(data_dir+tag+'_array.log',usecols = (2))
    lofar_trigg = np.loadtxt(data_dir+tag+'_array.log',usecols = (8))

    return timestamps, ns_timestamps, n_actived_detectors, lofar_trigg

def old_file_all_events_in_root(filename, data_dir, min_detectors):

    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_event = root_file.Get("Tree_event")
    num_entries = tree_event.GetEntries() 
    event_indices = [] 

    # Loop through all events
    for entry in range(num_entries):
        tree_event.GetEntry(entry) 
        
        active_detectors = 0 
        
        # Loop through all detectors
        for i in range(LORA.nDetA): 
            detname = f'Det{1 + i}'
            det_branch = tree_event.GetBranch(detname)
            nsec_leaf = det_branch.GetLeaf("nsec")
            nsec_value = nsec_leaf.GetValue()
                    
            if nsec_value > 0: 
                active_detectors += 1
                if active_detectors >= min_detectors:
                    event_indices.append(entry)
                    break 
    timestamps = []
    nstimestamps = []
    for e in range(len(event_indices)):
        entry = event_indices[e]
        tree_event.GetEntry(int(entry)) 
        
        gps_time, nsec_time = None, None
        
        # Loop through all detectors to find the first valid timestamp
        for i in range(LORA.nDetA): 
            detname = f'Det{1 + i}'
            det_branch = tree_event.GetBranch(detname)
            gps_leaf = det_branch.GetLeaf("GPS_time_stamp")
            nsec_leaf = det_branch.GetLeaf("nsec")
            gps_value = gps_leaf.GetValue()
            nsec_value = nsec_leaf.GetValue()
            if nsec_value > 0:
                gps_time, nsec_time = gps_value, nsec_value
                break
        
        
        timestamps.append(gps_time)
        nstimestamps.append(nsec_time)
    return timestamps, nstimestamps

def num_triggered_det(min_detectors, trigger_array):

    return event_indices


def load_in_data_V1(LOFAR_id, data_dir):
    #Generate timestamp
    timestamp=int(float(LOFAR_id))+LORA.event_id_offset
    print(timestamp)
    #Load in detectors
    detectors=[]
    lasas=[]

    for i in np.arange(LORA.nDetA):
        detectors.append(detector.Detector('Det'+(str(i+1))))

    for i in np.arange(LORA.nLasaA):
        lasas.append(detector.Lasa('Lasa'+(str(i+1))))

    #Load in position, signal and gain                           
    detector.load_positions(detectors)
    detector.load_signal(detectors)
    detector.load_gain(detectors)

    #Find timestamp in data files
    LORA4times=np.genfromtxt('/vol/astro5/lofar/astro3/vhecr/lora_triggered/LORA/'+'LORAtime4')
    index=-1

    for i in np.arange(len(LORA4times)):
        if int(LORA4times[i][0])== timestamp:
            index=i
            continue
    
    #Check if timestamp exists
    if i>-1:
        run=1
    else:
        run=0
        print('can\'t  find timestamp')

    #Load in event
    ns_timestamp = int(LORA4times[index][1])
    ev=event.Event(LOFAR_id,'V1')

    #Find tag of root file for easier data loading
    tag=find_tag(timestamp,ns_timestamp,data_dir)
    print(tag, timestamp)
    print('_____________________')
    print('{2}:    {0}   {1}'.format(timestamp,ns_timestamp,tag))

    #Use a longer method if the root does not have a tag
    if tag=='no_match':
        print('no log file found, doing this the hard way')
        tag=find_tag_exception(timestamp,ns_timestamp,data_dir)
        print('---> {0}'.format(tag))

    #Get general, log, time and noise data from root files
    log_file_info=log_file(tag,data_dir)
    info=return_root(tag,timestamp,ns_timestamp,data_dir)
    sec_info0,sec_info1,sec_info2=return_second_data(tag,timestamp,ns_timestamp,data_dir)
    log_info=return_log_data(tag,timestamp,ns_timestamp,data_dir)
    noise_info=return_noise_data(tag,timestamp,ns_timestamp,data_dir)

    #Load in trigger info
    LOFAR_trig=log_file_info['LOFAR_trig']
    lasa1_status=log_file_info['lasa1_status']
    lasa2_status=log_file_info['lasa2_status']
    lasa3_status=log_file_info['lasa3_status']
    lasa4_status=log_file_info['lasa4_status']
    lasa5_status=log_file_info['lasa5_status']

    

    #Load in event,sec, log and noise data into detector objects
    if (sec_info0[0]['GPS_time_stamp']-info[0]['gps'])>1:
        print('Can\'t find corresponding 1 Sec message: {0}  {1}'.format(int(info[0]['gps']),int(sec_info0[0]['GPS_time_stamp'])))
        for l in np.arange(5):
            lasas[l].sec_flag=1
    detector.load_event_information(info,detectors)
    detector.load_sec_information(sec_info0,sec_info1,sec_info2,lasas,'V1')
    detector.load_log_information(log_info,detectors)
    detector.load_noise_information(noise_info,detectors)

    #Get arrival times and counts in detectors 
    for d in np.arange(LORA.nDetA):
        #print('trigg_condition'+str(d)+':'+str(detectors[d].trigg_condition))
        event.find_counts(detectors[d])
        event.get_arrival_time(detectors[d])
        event.get_event_timestamp(detectors[d],lasas[int((detectors[d].number-1)/4)])
        event.find_counts_backgroundV2(detectors[d])
        detectors[d].trace_int_counts=detectors[d].total_counts_backgroundV2
        #print('threshold: {0:.2f},  corr peak: {1:.2f}, baseline: {2:.2f}, mean: {6:.2f}, sigma: {7:.2f}, peak: {3:.2f},  trigger: {4}, {5}'.format(detectors[d].threshold,detectors[d].corrected_peak,detectors[d].trace_mean, detectors[d].peak,detectors[d].trig,detectors[d].trigg_condition,detectors[d].sec_mean,detectors[d].sec_sigma))

    #Calculate event timestamp for further time calculations
    for l in np.arange(LORA.nLasaA):
        event.cal_event_timestamp(detectors,lasas[l])

    return detectors, ev

def load_in_data_V2(LOFAR_id, data_dir):
    
    #Generate timestamp
    timestamp=int(float(LOFAR_id))+LORA.event_id_offset

    #Load in detectors
    detectors=[]
    lasas=[]

    for i in np.arange(LORA.nDetB):
        detectors.append(detector.Detector('Det'+(str(i+1))))

    for i in np.arange(LORA.nLasaB):
        lasas.append(detector.Lasa('Lasa'+(str(i+1))))

    #Load in positions, signal and gains                           
    detector.load_positions(detectors)
    detector.load_signal(detectors)
    detector.load_gain(detectors)

    #Find timestamp in data files
    LORA4times=np.genfromtxt('/vol/astro5/lofar/astro3/vhecr/lora_triggered/LORA/'+'LORAtime4')
    index=-1

    for i in np.arange(len(LORA4times)):
        if int(LORA4times[i][0])== timestamp:
            index=i
            continue
    
    #Check if timestamp exists
    if i>-1:
        run=1
    else:
        run=0
        print('can\'t  find timestamp')

    #Load in event
    ns_timestamp = int(LORA4times[index][1])
    ev=event.Event(LOFAR_id,'V2')

    #Find tag of root file for easier data laoding
    tag=find_tag(timestamp,ns_timestamp,data_dir)
    print(tag, timestamp)
    print('_____________________')
    print('{2}:    {0}   {1}'.format(timestamp,ns_timestamp,tag))

    #Use a longer method if the root does not have a tag
    if tag=='no_match':
        print('no log file found, doing this the hard way')
        tag=find_tag_exceptionV2(timestamp,ns_timestamp,data_dir)
        print('---> {0}'.format(tag))

    #Get general, log, time and noise data from root files
    log_file_info=log_file(tag,data_dir)
    info=return_rootV2(tag,timestamp,ns_timestamp,data_dir)
    sec_info0,sec_info1,sec_info2=return_second_dataV2(tag,timestamp,ns_timestamp,data_dir)
    log_info=return_log_dataV2(tag,timestamp,ns_timestamp,data_dir)
    noise_info=return_noise_dataV2(tag,timestamp,ns_timestamp,data_dir)

    #Load in trigger info
    LOFAR_trig=log_file_info['LOFAR_trig']
    lasa1_status=log_file_info['lasa1_status']
    lasa2_status=log_file_info['lasa2_status']
    lasa3_status=log_file_info['lasa3_status']
    lasa4_status=log_file_info['lasa4_status']
    lasa5_status=log_file_info['lasa5_status']

    

    #Load in event,sec, log and noise data into detector objects
    if (sec_info0[0]['GPS_time_stamp']-info[0]['gps'])>1:
        #print('Can\'t find corresponding 1 Sec message: {0}  {1}'.format(int(info[0]['gps']),int(sec_info0[0]['GPS_time_stamp'])))
        for l in np.arange(10):
            lasas[l].sec_flag=1
    detector.load_event_information(info,detectors)
    detector.load_sec_informationV2(sec_info0,sec_info1,sec_info2,lasas,'V2')
    detector.load_log_information(log_info,detectors)
    detector.load_noise_information(noise_info,detectors)

    #Get arrival times and counts in detectors
    for d in np.arange(LORA.nDetB):
        event.find_counts(detectors[d])
        event.get_arrival_time(detectors[d])
        event.get_event_timestamp_V2(detectors[d],lasas[int((detectors[d].number-1)/4)])
        #event.find_counts_backgroundV2(detectors[d])
        #detectors[d].trace_int_counts=detectors[d].total_counts_backgroundV2
        #print('threshold: {0:.2f},  corr peak: {1:.2f}, baseline: {2:.2f}, mean: {6:.2f}, sigma: {7:.2f}, peak: {3:.2f},  trigger: {4}, {5}'.format(detectors[d].threshold,detectors[d].corrected_peak,detectors[d].trace_mean, detectors[d].peak,detectors[d].trig,detectors[d].trigg_condition,detectors[d].sec_mean,detectors[d].sec_sigma))

    #Calculate event timestamp for further time calculations
    for l in np.arange(LORA.nLasaB):
        event.cal_event_timestamp(detectors,lasas[l])

    return detectors, ev




def find_tag(timestamp,ns_timestamp,data_dir):
    #print '\n________________________\n'
    #timestamp = LORA4times[t][0]
    #ns_timestamp = LORA4times[t][1]
    
    dt_object = datetime.fromtimestamp(timestamp)
    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    second=dt_object.second
    
    #print year,month,day,hour
    #print int(timestamp),int(ns_timestamp)
   

    filestr1= str(year)+str(month).zfill(2)+str(day).zfill(2)
    filestr3='-1'
    filestr4='-1'
    
    if day>1:
        filestr2= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
        filestr4= str(year)+str(month).zfill(2)+str(day-2).zfill(2)
    
    if (month==10 or month==5 or month==7 or month==12) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
    
    if (month==2 or month==4 or month==6 or month==8 or month==9 or month==11) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(31).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
    
    if (month==3) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(28).zfill(2)
        filestr3= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(27).zfill(2)
    
    if (month==1) and day==1:
        filestr2= str(year-1)+str(12).zfill(2)+str(31).zfill(2)
        filestr4= str(year-1)+str(12).zfill(2)+str(30).zfill(2)

    #filestr3= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
    
    filelist=[]
    
    for file in os.listdir(data_dir):
        if filestr1 in file or filestr2 in file or filestr3 in file or filestr4 in file:
            if '.log' in file:
                #print(os.path.join(data_dir, file))
                filelist.append(os.path.join(data_dir, file))


    found=0
    filetag='no'
    for i in np.arange(len(filelist)):
        #print filelist[i]
        #for i in np.arange(1):
        file=open(filelist[i])
        stamps=[]
        try:
            stamps=np.genfromtxt(file,skip_header=10,usecols=(2,3))
        except:
            print('empty log')
        file.close()

        if len(stamps)>0:
            if len(stamps[(stamps.T[0]==timestamp)*(stamps.T[1]==ns_timestamp)])>0:
                
                filetag=filelist[i].split('raw/')[1].split('.')[0]
                found=1


    #print 'looking for this file: ',data_dir+filetag+'.root'

    if found==1 and os.path.exists(data_dir+filetag+'.root'):
        return filetag
    else:
        return 'no_match'

def find_tag_exception(timestamp,ns_timestamp,data_dir):
    #print '\n________________________\n'
    #timestamp = LORA4times[t][0]
    #ns_timestamp = LORA4times[t][1]
    
    print('no standard log file avaliable, trying to look at .root files')
    
    dt_object = datetime.fromtimestamp(timestamp)
    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    second=dt_object.second
    
    #print year,month,day,hour
    #print int(timestamp),int(ns_timestamp)
    filestr1= str(year)+str(month).zfill(2)+str(day).zfill(2)
    filestr3='-1'
    filestr4='-1'
    
    if day>1:
        filestr2= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
        filestr4= str(year)+str(month).zfill(2)+str(day-2).zfill(2)

    if (month==10 or month==5 or month==7 or month==12) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(29).zfill(2)

    if (month==2 or month==4 or month==6 or month==8 or month==9 or month==11) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(31).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(30).zfill(2)

    if (month==3) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(28).zfill(2)
        filestr3= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(27).zfill(2)

    if (month==1) and day==1:
        filestr2= str(year-1)+str(12).zfill(2)+str(31).zfill(2)
        filestr4= str(year-1)+str(12).zfill(2)+str(30).zfill(2)




    filelist=[]
    
    for file in os.listdir(data_dir):
        if filestr1 in file or filestr2 in file or filestr3 in file or filestr4 in file:
            if '.root' in file:
                #print(os.path.join(data_dir, file))
                filelist.append(os.path.join(data_dir, file))
    
    found=0
    filetag='no'

    for i in np.arange(len(filelist)):
        root_file=ROOT.TFile.Open(filelist[i])
        try:
            tree_event = root_file.Get("Tree_event")
        
        
            event_index=-1
            event_index=find_entry_number(timestamp,ns_timestamp,tree_event)
            print('did we find the index?   {0}'.format(event_index))
        
            if event_index>-1:
                found=1
                filetag=filelist[i].split('raw/')[1].split('.')[0]
                break
        except:
            print('can\'t open the event leaf')
    
    if found==1:
        print('returning file tag')
        return filetag
    
    else:
        return 'no_match'

def find_tag_exceptionV2(timestamp,ns_timestamp,data_dir):
    #print '\n________________________\n'
    #timestamp = LORA4times[t][0]
    #ns_timestamp = LORA4times[t][1]
    
    print('no standard log file avaliable, trying to look at .root files')
    
    dt_object = datetime.fromtimestamp(timestamp)
    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    second=dt_object.second
    
    #print year,month,day,hour
    #print int(timestamp),int(ns_timestamp)
    filestr1= str(year)+str(month).zfill(2)+str(day).zfill(2)
    filestr3='-1'
    filestr4='-1'
    
    if day>1:
        filestr2= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
        filestr4= str(year)+str(month).zfill(2)+str(day-2).zfill(2)

    if (month==10 or month==5 or month==7 or month==12) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(29).zfill(2)

    if (month==2 or month==4 or month==6 or month==8 or month==9 or month==11) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(31).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(30).zfill(2)

    if (month==3) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(28).zfill(2)
        filestr3= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(27).zfill(2)

    if (month==1) and day==1:
        filestr2= str(year-1)+str(12).zfill(2)+str(31).zfill(2)
        filestr4= str(year-1)+str(12).zfill(2)+str(30).zfill(2)




    filelist=[]
    
    for file in os.listdir(data_dir):
        if filestr1 in file or filestr2 in file or filestr3 in file or filestr4 in file:
            if '.root' in file:
                #print(os.path.join(data_dir, file))
                filelist.append(os.path.join(data_dir, file))
    
    found=0
    filetag='no'
    for i in np.arange(len(filelist)):
        root_file=ROOT.TFile.Open(filelist[i])
        #try:
        tree_event = root_file.Get("Tree_Event_Header")
    
    
        event_index=-1
        event_index=find_entry_numberV2(timestamp,ns_timestamp,tree_event)
        print('did we find the index?   {0}'.format(event_index))
    
        if event_index>-1:
            found=1
            filetag=filelist[i].split('raw/')[1].split('.')[0]
            break
        #except:
            print('can\'t open the event leaf')
    
    if found==1:
        print('returning file tag')
        return filetag
    
    else:
        return 'no_match'

    
def getTime(det, entry):
    det.GetEntry(entry)
    ymd=det.GetLeaf('YMD').GetValue()
    gps=det.GetLeaf('GPS_time_stamp').GetValue()
    ctd=det.GetLeaf('CTD').GetValue()
    nsec=det.GetLeaf('nsec').GetValue()
    return ymd,gps,ctd,nsec

def getTimeSec(lasa, entry):
    lasa.GetEntry(int(entry))
    gps=lasa.GetLeaf('GPS_time_stamp').GetValue()
    return gps

def getTimeSecV2(lasa, entry):
    lasa.GetEntry(int(entry))
    gps=lasa.GetLeaf('GPS_Time_Stamp').GetValue()
    return gps

def find_entry_number(lora_utc,lora_nsec,tree_event):
    event=-1
    
    det1=tree_event.GetBranch('Det1')
    det5=tree_event.GetBranch('Det5')
    det9=tree_event.GetBranch('Det9')
    det13=tree_event.GetBranch('Det13')
    det17=tree_event.GetBranch('Det17')
    det1.GetLeaf('GPS_time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
    trigger_check=0
    diff_best=1e10
    for e in range(nE):
        ymd1,gps1,ctd1,nsec1= getTime(det1,e)
        ymd2,gps2,ctd2,nsec2= getTime(det5,e)
        ymd3,gps3,ctd3,nsec3= getTime(det9,e)
        ymd4,gps4,ctd4,nsec4= getTime(det13,e)
        ymd5,gps5,ctd5,nsec5= getTime(det17,e)
        times=[gps1,gps2,gps3,gps4,gps5]
        times_ns=[nsec1,nsec2,nsec3,nsec4,nsec5]
        
        if lora_utc in times:
            diff=np.max(np.abs(lora_nsec-np.asarray(times_ns)[np.asarray(times_ns)>1]))
            if diff<10000:  # w/in 10 us to avoid mis-triggers
                if diff<diff_best:
                    diff_best=diff
                    trigger_check=1
                    event=e



    return event

def find_entry_numberV2(lora_utc,lora_nsec,tree_event):
    event=-1
    GPS=tree_event.GetLeaf('GPS_Time_Stamp_FirstHit')
    ns=tree_event.GetLeaf('nsec_Online_FirstHit')
    nE= tree_event.GetEntries()
    #times=np.zeros([nLasa,nE])
    trigger_check=0
    diff_best=1e10

    
    for e in range(nE):
        tree_event.GetEntry(e)
        times = GPS.GetValue() 
        times_ns = ns.GetValue()
        if lora_utc == times:
            diff=np.abs(lora_nsec-times_ns)
            if diff<10000:  # w/in 10 us to avoid mis-triggers
                if diff<diff_best:
                    diff_best=diff
                    trigger_check=1
                    event=e



    return event



def getDataV1(det, entry):
    
    det.GetEntry(int(entry))
    
    detector=det.GetLeaf('detector').GetValue()
    ymd=det.GetLeaf('YMD').GetValue()
    gps=det.GetLeaf('GPS_time_stamp').GetValue()
    ctd=det.GetLeaf('CTD').GetValue()
    nsec=det.GetLeaf('nsec').GetValue()
    trigg_condition=det.GetLeaf('Trigg_condition').GetValue()
    try:
        trigg_pattern=det.GetLeaf('Trigg_pattern').GetValue()
    except:
        trigg_pattern=-1
    total_counts=det.GetLeaf('Total_counts').GetValue()
    pulse_height=det.GetLeaf('Pulse_height').GetValue()
    pulse_width=det.GetLeaf('Pulse_width').GetValue()
    counts=det.GetLeaf('counts')
    hold=np.zeros([nTraceV1])
    for i in range(nTraceV1):
        hold[i]=counts.GetValue(i)

    info={'det':detector,'ymd':ymd,'gps':gps,'ctd':ctd,'nsec':nsec,'trigg_condition':trigg_condition,'trigg_pattern':trigg_pattern,'total_counts':total_counts,'pulse_height':pulse_height,'pulse_width':pulse_width,'counts':hold}
    return info

def getDataV2(tree, entry):
    if entry >= 0:
        tree.GetEntry(int(entry))
        
        detector=tree.GetLeaf('Detector').GetValue()
        #ymd=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
        ymd = 0
        gps=tree.GetLeaf('GPS_Time_Stamp').GetValue()
        ctd=tree.GetLeaf('CTD').GetValue()
        nsec=tree.GetLeaf('nsec_Online').GetValue()
        trigg_condition=tree.GetLeaf('HiSparc_Trigg_Condition').GetValue()
        try:
            trigg_pattern=tree.GetLeaf('HiSparc_Trigg_Pattern').GetValue()
        except:
            trigg_pattern=-1
        total_counts=tree.GetLeaf('Charge_Corrected').GetValue()
        pulse_height=tree.GetLeaf('Peak_Height_Corrected').GetValue()
        #pulse_width=tree.GetLeaf('Pulse_width').GetValue()
        counts=tree.GetLeaf('Waveform_Raw')
        hold=np.zeros([nTraceV1])
        for i in range(nTraceV1):
            hold[i]=counts.GetValue(i)

        info={'det':detector,'ymd':ymd,'gps':gps,'ctd':ctd,'nsec':nsec,'trigg_condition':trigg_condition,'trigg_pattern':trigg_pattern,'total_counts':total_counts,'pulse_height':pulse_height,'pulse_width':0,'counts':hold}
    
    else:
        info = {'det':0,'ymd':0,'gps':0,'ctd':0,'nsec':0,'trigg_condition':0,'trigg_pattern':0,'total_counts':0,'pulse_height':0,'pulse_width':0,'counts':np.zeros([nTraceV1])}
    return info


def get_sample_files(begin_date, data_dir, time_frame):

    start_dt = datetime.strptime(begin_date, "%Y%m%d")

    
    # Generate a set of allowed date prefixes (including the next 50 days)
    allowed_dates = { (start_dt + timedelta(days=i)).strftime("%Y%m%d") for i in range(time_frame+1) }
    
    files = os.listdir(data_dir)

    matching_files = sorted([f for f in files if any(f.startswith(date) for date in allowed_dates) and f.endswith('.root')])
    
    return matching_files

def sample_detector_output_average(filenames,data_dir, detector_number):
    nCounts = 400
    total_times = []
    avg_rms = []
    avg_mean = []
    
    for filename in filenames:
        try:
            root_file = ROOT.TFile.Open(data_dir + filename)
            tree_ev = root_file.Get("Tree_Event")
            n_entries = tree_ev.GetEntries()
            
            total_rms = []
            total_means = []

            for e in range(n_entries):
                tree_ev.GetEntry(e)
                detector_id = tree_ev.GetLeaf("Detector").GetValue()
                if detector_id == detector_number:
                    
                    counts_leaf = tree_ev.GetLeaf('Waveform_Raw')

                    counts = np.array([counts_leaf.GetValue(i) for i in range(nCounts)])
                    mean = np.mean(counts)
                    rms = np.sqrt(np.mean(np.square(counts)))

                    total_rms.append(rms)
                    total_means.append(mean)
            
            date = filename[0:8]
            total_times.append(date)
        
            avg_rms.append(np.mean(total_rms))
            avg_mean.append(np.mean(total_means))

        except OSError as error :
            print('root file not found, moving on') 
        except AttributeError as e:
            print(f"AttributeError occurred: {e}")
    total_times_formatted = [datetime.strptime(d, "%Y%m%d").strftime("%Y/%m/%d") for d in total_times]
    return avg_mean, avg_rms, total_times_formatted


def sample_detector_output_per_day(filenames,data_dir, detector_number):
    nCounts = 400
    total_times = []
    total_rms = []
    total_means = []
    
    for filename in filenames:
        root_file = ROOT.TFile.Open(data_dir + filename)
        tree_ev = root_file.Get("Tree_Event")
        n_entries = tree_ev.GetEntries()
        found_entry = False 

        for e in range(n_entries):
            tree_ev.GetEntry(e)
            detector_id = tree_ev.GetLeaf("Detector").GetValue()
            if detector_id == detector_number:
                date = filename[0:8]
                counts_leaf = tree_ev.GetLeaf('Waveform_Raw')

                counts = np.array([counts_leaf.GetValue(i) for i in range(nCounts)])
                mean = np.mean(counts)
                rms = np.sqrt(np.mean(np.square(counts)))

                total_times.append(date)
                total_rms.append(rms)
                total_means.append(mean)

                found_entry = True
                break  

        if not found_entry:
            print(f"Warning: No valid entry found for detector {detector_number} in {filename}")

    total_times_formatted = [datetime.strptime(d, "%Y%m%d").strftime("%Y/%m/%d") for d in total_times]
    return total_means, total_rms, total_times_formatted

def GetRMS(source_dir, file_name):
    file_path = source_dir + file_name + '.dat'

    # Read file line by line
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip header

    detector_data = {}

    for line in lines:
        cols = line.strip().split("\t")  # Assuming tab-separated values

        if len(cols) < 4:
            continue  # Skip malformed lines

        try:
            det = int(cols[0])  # Detector ID
            rms = float(cols[1]) if cols[1] != "nan" else np.nan
            mean = float(cols[2]) if cols[2] != "nan" else np.nan
            date = cols[3]  # Keep date as a string

            # Store in dictionary
            if det not in detector_data:
                detector_data[det] = ([], [], [])

            detector_data[det][0].append(rms)   # Append RMS
            detector_data[det][1].append(mean)  # Append Mean
            detector_data[det][2].append(date)  # Append Date

        except ValueError:
            print(f"Skipping malformed line: {line}")
            continue  # Skip problematic lines

    return detector_data

def find_data_numberV2(lora_utc,tree_noise,d):
    
    det1=tree_noise.GetLeaf('UTC_Time_stamp')
    nE= tree_noise.GetEntries()
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        tree_noise.GetEntry(int(e))
        gps1=tree_noise.GetLeaf('GPS_Time_Stamp').GetValue()
        detector_id = tree_noise.GetLeaf('Detector').GetValue()
        # Match the detector ID and the timestamp
        if gps1 == lora_utc and detector_id == d and index_found == 0:
            event = e
            index_found = 1
    return event

def return_root(filename,utc,nsec,data_dir):

    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    print('reading root file: {0}'.format(data_dir+filename+'.root'))
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")
    event_index=find_entry_number(utc,nsec,tree_event)
    all_info=[]
    for i in np.arange(LORA.nDetA):
        detname='Det'+str(1+i)
        det=tree_event.GetBranch(detname)   
        info=getDataV1(det,event_index)
        all_info.append(info)
    

    return all_info


def return_rootV2(filename,utc,nsec,data_dir):

    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    print('reading root file: {0}'.format(data_dir+filename+'.root'))
    tree_sec = root_file.Get("Tree_sec")
    tree_ev=root_file.Get("Tree_Event")
    tree_event = root_file.Get("Tree_Event_Header")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")
    event_index=find_entry_numberV2(utc,nsec,tree_event)
    entries = []
    for i in np.arange(LORA.nDetB):
        d = i+1
        entries.append(find_data_numberV2(utc, tree_ev, d))


    all_info=[]
    for i in np.arange(LORA.nDetB):
        info=getDataV2(tree_ev, entries[i])
        all_info.append(info)
    
    return all_info

def find_sec_number(lora_utc,tree_sec,i):
    event=-1
    
    lasa1=tree_sec.GetBranch('Lasa'+str(i+1))


    lasa1.GetLeaf('GPS_time_stamp')
    nE= lasa1.GetEntries()
    #times=np.zeros([nLasa,nE])
    diff_best=1e10
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        gps1= getTimeSec(lasa1,e)
        if gps1==lora_utc and index_found==0:
            event=e
            index_found=1
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1

    return event


def find_sec_numberV2(lora_utc,tree_sec,i):
    event=-1
    
    lasa1 = tree_sec.GetLeaf('GPS_Time_Stamp')
    nE = tree_sec.GetEntries()
    #times=np.zeros([nLasa,nE])
    diff_best=1e10
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        gps1= getTimeSecV2(tree_sec,e)
        if gps1==lora_utc and index_found==0:
            event=e
            index_found=1
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1

    return event


def getSecV1(det, entry):
    
    det.GetEntry(entry)
    
    lasa=det.GetLeaf('Lasa').GetValue()
    
    YMD=det.GetLeaf('YMD').GetValue()

    GPS_time_stamp=det.GetLeaf('GPS_time_stamp').GetValue()
    sync=det.GetLeaf('sync').GetValue()
    CTP=det.GetLeaf('CTP').GetValue()
    quant=det.GetLeaf('quant').GetValue()
    Channel_1_Thres_count_high=det.GetLeaf('Channel_1_Thres_count_high').GetValue()
    Channel_1_Thres_count_low=det.GetLeaf('Channel_1_Thres_count_low').GetValue()
    Channel_2_Thres_count_high=det.GetLeaf('Channel_2_Thres_count_high').GetValue()
    Channel_2_Thres_count_low=det.GetLeaf('Channel_2_Thres_count_low').GetValue()
    Satellite_info=det.GetLeaf('Satellite_info').GetValue()

    info={'lasa':lasa,'YMD':YMD,'GPS_time_stamp':GPS_time_stamp,'sync':sync,'CTP':CTP,'quant':quant}
    return info

def getSecV2(det, entry):
    
    det.GetEntry(entry)
    
    lasa=det.GetLeaf('Station').GetValue()
    
    #YMD=det.GetLeaf('YMD').GetValue()

    GPS_time_stamp=det.GetLeaf('GPS_Time_Stamp').GetValue()
    sync=det.GetLeaf('Sync_Error').GetValue()
    CTP=det.GetLeaf('CTP').GetValue()
    quant=det.GetLeaf('Quant_Error').GetValue()
    Channel_1_Thres_count_high=det.GetLeaf('Channel_1_Thres_Count_High').GetValue()
    Channel_1_Thres_count_low=det.GetLeaf('Channel_1_Thres_Count_Low').GetValue()
    Channel_2_Thres_count_high=det.GetLeaf('Channel_2_Thres_Count_High').GetValue()
    Channel_2_Thres_count_low=det.GetLeaf('Channel_2_Thres_Count_Low').GetValue()
    Satellite_info=det.GetLeaf('Satellite_Info').GetValue()

    info={'lasa':lasa,'YMD':0,'GPS_time_stamp':GPS_time_stamp,'sync':sync,'CTP':CTP,'quant':quant}
    return info


def return_second_data(filename,utc,nsec,data_dir):

    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")

    entry=np.zeros([LORA.nLasaA])
    
    for i in np.arange(LORA.nLasaA):
        entry[i]=find_sec_number(utc,tree_sec,i)
    
    all_info=[]
    all_info1=[]
    all_info2=[]

    
    for i in np.arange(LORA.nLasaA):
        lasaname='Lasa'+str(1+i)
        det=tree_sec.GetBranch(lasaname)
        all_info.append(getSecV1(det, int(entry[i])))
        all_info1.append(getSecV1(det, int(entry[i]+1)))
        all_info2.append(getSecV1(det, int(entry[i]+2)))

    return all_info,all_info1,all_info2

def return_second_dataV2(filename,utc,nsec,data_dir):

    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_Event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")

    entry=np.zeros([LORA.nLasaB])
    
    for i in np.arange(LORA.nLasaB):
        entry[i]=find_sec_numberV2(utc,tree_event,i)
    
    all_info=[]
    all_info1=[]
    all_info2=[]

    
    for i in np.arange(LORA.nLasaB):
        lasaname='Lasa'+str(1+i)
        dat=root_file.Get('Tree_OSM_HiSparc')
        all_info.append(getSecV2(dat, int(entry[i])))
        all_info1.append(getSecV2(dat, int(entry[i]+1)))
        all_info2.append(getSecV2(dat, int(entry[i]+2)))

    return all_info,all_info1,all_info2


def getLogV1(det,d, entry):
    nE= det.GetEntries()
    det.GetEntry(int(entry))
    
    
    YMD=det.GetLeaf('YMD').GetValue()
    GPS_time_stamp=det.GetLeaf('Time_stamp').GetValue()
    Threshold_low=det.GetLeaf('Channel_thres_low').GetValue()
    
    #print('_________________________________________')
    #print('finding threshold:  {0}'.format(Threshold_low))
    if Threshold_low==0.0:
        print('issue with Threshold_low {0}'.format(nE))
        thresh_avg=0
        thresh_count=0
        for i in np.arange(nE):
            det.GetEntry(int(i))
            thesh_temp=det.GetLeaf('Channel_thres_low').GetValue()
            if thesh_temp>0.0:
                thresh_avg=thresh_avg+thesh_temp
                thresh_count=thresh_count+1
        if thresh_count>0:
            Threshold_low=thresh_avg/(1.0*thresh_count)
        else:
            Threshold_low=avg_thresholdV1[d]
        
    
    info={'threshold':Threshold_low}
    return info

def getLogV2(det, d, entry):
    nE = det.GetEntries()
    det.GetEntry(int(entry))

    Threshold_low = det.GetLeaf('Channel_Thres_Low').GetValue()

    #print('_________________________________________')
    #print(f'Finding threshold for detector {d}: {Threshold_low}')

    # Handle cases where the threshold is zero
    if Threshold_low == 0.0:
        print(f'Issue with Threshold_low for detector {d}, calculating average.')
        thresh_avg = 0
        thresh_count = 0
        for i in np.arange(nE):
            det.GetEntry(i)
            detector_id_temp = det.GetLeaf('Detector').GetValue()
            if detector_id_temp == d:
                thesh_temp = det.GetLeaf('Channel_Thres_Low').GetValue()
                if thesh_temp > 0.0:
                    thresh_avg += thesh_temp
                    thresh_count += 1
        if thresh_count > 0:
            Threshold_low = thresh_avg / thresh_count
        else:
            Threshold_low = avg_thresholdV1[d]  # Fallback to predefined average if necessary

    return {'threshold': Threshold_low}


def getNoiseV1(det, entry):
    
    det.GetEntry(int(entry))
    
    sigma=det.GetLeaf('Sigma').GetValue()
    mean=det.GetLeaf('Mean').GetValue()
    info={'mean':mean,'sigma':sigma}
    return info

def getNoiseV2(det, entry):
    
    det.GetEntry(int(entry))
    
    sigma=det.GetLeaf('Mean_Sigma').GetValue()
    mean=det.GetLeaf('Mean_Baseline').GetValue()
    info={'mean':mean,'sigma':sigma}
    return info


def find_noise_number(lora_utc,tree_noise,d):
    
    event=-1
    
    det1=tree_noise.GetBranch('Det'+str(d))
    
    det1.GetLeaf('GPS_time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
    
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        det1.GetEntry(int(e))
        gps1=det1.GetLeaf('GPS_time_stamp').GetValue()
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1
    return event

def find_noise_numberV2(lora_utc,tree_noise,d):
    
    det1=tree_noise.GetLeaf('UTC_Time_stamp')
    nE= tree_noise.GetEntries()
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        tree_noise.GetEntry(int(e))
        gps1=tree_noise.GetLeaf('UTC_Time_Stamp').GetValue()
        detector_id = tree_noise.GetLeaf('Detector').GetValue()
        
        # Match the detector ID and the timestamp
        if gps1 > lora_utc and detector_id == d and index_found == 0:
            event = e
            index_found = 1
    return event

def find_log_number(lora_utc,tree_log,d):
    
    event=-1
    
    det1=tree_log.GetBranch('Det'+str(d))
    
    det1.GetLeaf('Time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
  
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        det1.GetEntry(int(e))
        gps1=det1.GetLeaf('Time_stamp').GetValue()
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1

    return event


def find_log_numberV2(lora_utc, tree_log, d):
    nE = tree_log.GetEntries()

    event = -1
    index_found = 0

    # Search for the correct entry in the log
    for e in np.arange(nE):
        tree_log.GetEntry(int(e))
        gps1 = tree_log.GetLeaf('UTC_Time_Stamp').GetValue()
        detector_id = tree_log.GetLeaf('Detector').GetValue()
        if gps1 > lora_utc and detector_id == d and index_found == 0:
            event = e
            index_found = 1

    return event

def return_log_data(filename,utc,nsec,data_dir):
    
    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")
    
    entry1=find_log_number(utc,tree_log,1)
    entry2=find_log_number(utc,tree_log,5)
    entry3=find_log_number(utc,tree_log,9)
    entry4=find_log_number(utc,tree_log,13)
    entry5=find_log_number(utc,tree_log,17)
    
    all_info=[]
    
    
    for i in np.arange(LORA.nDetB):
        detname='Det'+str(1+i)
        det=tree_log.GetBranch(detname)
        if i==0 or i==1 or i==2 or i==3:
            all_info.append(getLogV1(det,i, entry1))
        if i==4 or i==5 or i==6 or i==7:
            all_info.append(getLogV1(det,i, entry2))
        if i==8 or i==9 or i==10 or i==11:
            all_info.append(getLogV1(det,i, entry3))
        if i==12 or i==13 or i==14 or i==15:
            all_info.append(getLogV1(det,i, entry4))
        if i==16 or i==17 or i==18 or i==19:
            all_info.append(getLogV1(det,i, entry5))
    
    return all_info

def return_log_dataV2(filename, utc, nsec, data_dir):
    # Open ROOT file and retrieve trees
    root_file = ROOT.TFile.Open(data_dir + filename + '.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_Log")
    tree_noise = root_file.Get("Tree_noise")

    # Find entries in the log tree corresponding to the UTC time
    entries = []
    for i in np.arange(LORA.nDetB):
        d = i+1
        entries.append(find_log_numberV2(utc, tree_log, d))
    all_info = []

    # Retrieve threshold information for each detector
    for i in np.arange(LORA.nDetB):
        d = i+1
        all_info.append(getLogV2(tree_log, d, entries[i]))

    return all_info





def return_noise_data(filename,utc,nsec,data_dir):
    
    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")
    
    
    entry1=int(find_noise_number(utc,tree_noise,1))
    entry2=int(find_noise_number(utc,tree_noise,5))
    entry3=int(find_noise_number(utc,tree_noise,9))
    entry4=int(find_noise_number(utc,tree_noise,13))
    entry5=int(find_noise_number(utc,tree_noise,17))
    all_info=[]
    
   
    for i in np.arange(LORA.nDetB):
        detname='Det'+str(1+i)
        det=tree_noise.GetBranch(detname)
        if i==0 or i==1 or i==2 or i==3:
            all_info.append(getNoiseV1(det, int(entry1)))
        if i==4 or i==5 or i==6 or i==7:
            all_info.append(getNoiseV1(det, int(entry2)))
        if i==8 or i==9 or i==10 or i==11:
            all_info.append(getNoiseV1(det, int(entry3)))
        if i==12 or i==13 or i==14 or i==15:
            all_info.append(getNoiseV1(det, int(entry4)))
        if i==16 or i==17 or i==18 or i==19:
            all_info.append(getNoiseV1(det, int(entry5)))
    
    return all_info

def return_noise_dataV2(filename,utc,nsec,data_dir):
    
    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_Log")
    tree_noise = root_file.Get("Tree_Noise")


    entries = []
    for i in np.arange(LORA.nDetB):
        d = i+1
        entries.append(find_noise_numberV2(utc, tree_log, d))

    all_info=[]
    
    for i in np.arange(LORA.nDetB):
        all_info.append(getNoiseV2(tree_log, entries[i]))
    return all_info

def log_file(filename,data_dir):
    filepath=data_dir+filename+'.log'
    #log_file=open(data_dir+filename+'.log','r')
    lasa1_status=-1
    lasa2_status=-1
    lasa3_status=-1
    lasa4_status=-1
    lasa5_status=-1
    LOFAR_trig=-1
    try:
        with open(filepath,'r') as fp:
            line = fp.readline()
            cnt = 1
            while line:
            #print("Line {}: {}".format(cnt, line.strip()))
                if 'LOFAR trigger settings' in line:
                    LOFAR_trig=int(line.strip().split(':')[1])
            
                if 'CS003:' in line:
                    lasa1_status=int(line.strip().split(':')[1])
                if 'CS004:' in line:
                        lasa2_status=int(line.strip().split(':')[1])
                if 'CS005:' in line:
                            lasa3_status=int(line.strip().split(':')[1])
                if 'CS006:' in line:
                    lasa4_status=int(line.strip().split(':')[1])
                if 'CS007:' in line:
                    lasa5_status=int(line.strip().split(':')[1])
            
                line = fp.readline()
                cnt += 1

    except:
        print('can\'t find log file')
    '''
    print 'LOFAR trigger: ',LOFAR_trig
    print 'lasa 1: ',lasa1_status
    print 'lasa 2: ',lasa2_status
    print 'lasa 3: ',lasa3_status
    print 'lasa 4: ',lasa4_status
    print 'lasa 5: ',lasa5_status
    '''
    info={'LOFAR_trig':LOFAR_trig,'lasa1_status':lasa1_status,'lasa2_status':lasa2_status,'lasa3_status':lasa3_status,'lasa4_status':lasa4_status,'lasa5_status':lasa5_status}
    return info






def return_event_V2(event_id,event_GPS, event_ns,event_data):

    Station=event_data['Station']
    Detector=event_data['Detector']
    Channel_Passed_Threshold=event_data['Channel_Passed_Threshold']
    Trigg_Threshold=event_data['Trigg_Threshold']
    Charge_Corrected=event_data['Charge_Corrected']
    Peak_Height_Corrected=event_data['Peak_Height_Corrected']
    Peak_Height_Raw=event_data['Peak_Height_Raw']
    Waveform_Raw=event_data['Waveform_Raw']
    Event_Id=event_data['Event_Id']
    Run_Id=event_data['Run_Id']
    GPS_Time_Stamp=event_data['GPS_Time_Stamp']
    CTD=event_data['CTD']
    nsec_Online=event_data['nsec_Online']
    HiSparc_Trigg_Pattern=event_data['HiSparc_Trigg_Pattern']
    HiSparc_Trigg_Condition=event_data['HiSparc_Trigg_Condition']


    dets=Detector[Event_Id==event_id]
    counts=Waveform_Raw[Event_Id==event_id]
    pulse_height=Peak_Height_Corrected[Event_Id==event_id]
    total_counts=Charge_Corrected[Event_Id==event_id]
    trigger_pattern=HiSparc_Trigg_Pattern[Event_Id==event_id]
    trigger_condition=HiSparc_Trigg_Condition[Event_Id==event_id]
    nsec=nsec_Online[Event_Id==event_id]
    ctd=CTD[Event_Id==event_id]
    gps=GPS_Time_Stamp[Event_Id==event_id]
    thresh=Trigg_Threshold[Event_Id==event_id]
    #print('threshold--> ',thresh)
    #print counts.shape
    all_info=[]
    log_all_info=[]
    

    for d in np.arange(40):
        if (d+1) in dets:
            ind=np.where(dets==(d+1))[0][0]
            #print d+1, dets[ind], total_counts[ind]
            timestamp = date.fromtimestamp(gps[ind])
            ymd=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
            info={'det':dets[ind],'ymd':ymd,'gps':gps[ind],'ctd':ctd[ind],'nsec':nsec[ind],'trigg_condition':trigger_condition[ind],'trigg_pattern':trigger_pattern[ind],'total_counts':total_counts[ind],'pulse_height':pulse_height[ind],'pulse_width':0,'counts':counts[ind]}
            log_info={'threshold':thresh[ind]}
            all_info.append(info)
            log_all_info.append(log_info)
        else:
            info={'det':d+1,'ymd':0,'gps':0,'ctd':0,'nsec':0,'trigg_condition':0,'trigg_pattern':0,'total_counts':0,'pulse_height':0,'pulse_width':0,'counts':np.zeros([4000])}
            log_info={'threshold':0}

            all_info.append(info)
            log_all_info.append(log_info)



    return all_info,log_all_info



def return_second_data_V2(event_id,event_GPS, event_ns,osm_data_hisparc,osm_data_aera):


    Station_H=osm_data_hisparc['Station']
    Master_Or_Slave_H=osm_data_hisparc['Master_Or_Slave']
    GPS_Time_Stamp_H=osm_data_hisparc['GPS_Time_Stamp']
    Sync_Error_H=osm_data_hisparc['Sync_Error']
    Quant_Error_H=osm_data_hisparc['Quant_Error']
    CTP_H=osm_data_hisparc['CTP']
    
    
    Station_A=osm_data_aera['Station']
    GPS_Time_Stamp_A=osm_data_aera['GPS_Time_Stamp']
    Sync_Error_A=osm_data_aera['Sync_Error']
    Quant_Error_A=osm_data_aera['Quant_Error']
    CTP_A=osm_data_aera['CTP']
    UTC_offset_A=osm_data_aera['UTC_offset']
    
    
    all_info=[]
    all_info1=[]
    all_info2=[]

    for t in np.arange(3):
        #for i in np.arange(1):
        for i in np.arange(LORA.nLasaB):
            lasa=i+1
            #master
            #print('getting OSM for {0}'.format(i+1))
            if lasa<=5:
                gpsM= GPS_Time_Stamp_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==0)]
                if len(gpsM>0):
                    syncM= Sync_Error_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==0)]
                    quantM= Quant_Error_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==0)]
                    ctpM= CTP_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==0)]
                    print(gpsM)
                    timestamp = date.fromtimestamp(gpsM[0])
                    ymdM=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
                    #slave
                    gpsS= GPS_Time_Stamp_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==1)]
                    syncS= Sync_Error_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==1)]
                    #print((GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1))
                    #print((Station_H==i+1))
                    #print((Master_Or_Slave_H==1))
                    #print(Master_Or_Slave_H)
                    quantS= Quant_Error_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==1)]
                    ctpS= CTP_H[(GPS_Time_Stamp_H==(event_GPS+t))*(Station_H==i+1)*(Master_Or_Slave_H==1)]
                    timestamp = date.fromtimestamp(gpsM[0])
                    ymdS=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
                
                    info={'lasa':lasa,'YMD_M':ymdM,'GPS_time_stamp_M':gpsM,'sync_M':syncM,'CTP_M':ctpM,'quant_M':quantM,'YMD_S':ymdS,'GPS_time_stamp_S':gpsS,'sync_S':syncS,'CTP_S':ctpS,'quant_S':quantS}
                else:
                    info={'lasa':lasa,'YMD_M':np.asarray([0]),'GPS_time_stamp_M':np.asarray([0]),'sync_M':np.asarray([0]),'CTP_M':np.asarray([0]),'quant_M':np.asarray([0]),'YMD_S':np.asarray([0]),'GPS_time_stamp_S':np.asarray([0]),'sync_S':np.asarray([0]),'CTP_S':np.asarray([0]),'quant_S':np.asarray([0])}
            
            else:
                
                gps= GPS_Time_Stamp_A[(GPS_Time_Stamp_A==(event_GPS+t))*(Station_A==i+1)]
                #print(i+1,t, gps)
                if len(gps>0):
                    sync= Sync_Error_A[(GPS_Time_Stamp_A==(event_GPS+t))*(Station_A==i+1)]
                    quant= Quant_Error_A[(GPS_Time_Stamp_A==(event_GPS+t))*(Station_A==i+1)]
                    ctp= CTP_A[(GPS_Time_Stamp_A==(event_GPS+t))*(Station_A==i+1)]
                    timestamp = date.fromtimestamp(gps[0])
                    
                    ymd=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
                    info={'lasa':lasa,'YMD':ymd,'GPS_time_stamp':gps,'sync':sync,'CTP':ctp,'quant':quant}
                else:
                    info={'lasa':lasa,'YMD':np.asarray([0]),'GPS_time_stamp':np.asarray([0]),'sync':np.asarray([0]),'CTP':np.asarray([0]),'quant':np.asarray([0])}
                
            
            if t==0:
                all_info.append(info)
            if t==1:
                all_info1.append(info)
            if t==2:
                all_info2.append(info)
    
    return all_info,all_info1,all_info2
