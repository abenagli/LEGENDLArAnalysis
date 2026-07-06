// ________________________________________ //
// Author: Henrique Souza
// Filename: adc_read_all_data.C
// Created: 2021
// ________________________________________ //
#include "/home/carla/Documents/ADC_Henrique/cold_box_analysis/class/MYCODES.h"


void adc_read_all_data(string datadir = "./", float intTime = 1000){
    
    system("rm files.log");
    int datadirlength = datadir.length();
    if(datadir[datadirlength-1] != '/'){
        datadir = datadir + "/";
    }
    string packdata = Form("ls -1 -v %s*.dat | sed -e 's/.dat$//' > files.log", datadir.c_str());
    // string packadata = Form("ls -1 %s*.dat | grep wave | sed -e 's/.dat$//' > files.log", datadir.c_str()); // to only check one channel
    // string packadata = Form("ls -1 %s*.dat | grep wave[2,3] | sed -e 's/.dat$//' > files.log", datadir.c_str()); // to only check channels 2 and 3
    system(packdata.c_str());
    
    Read r;
        
    r.dtime = 4; // step time in ns
    r.nbits = 14; // this is used only for baseline histogram, changing it to 12 might help
    r.isBinary = true;
    
    r.baselineTime = 4000; // time limit for baseline in ns
    r.startCharge = 4000;  // start integration time in ns
    r.chargeTime = r.startCharge + intTime; // last time to integrate ns
    r.maxRange = 1000; // max range to search for amplitude peak
    r.fast = 600; // fprompt fast integration time
    r.slow = 3500; //fprompt slow integration time
    r.exclusion_baseline = 5; // filtered waveform (ADCs), anything above here will do +exclusion window while evaluating the baseline
    // r.exclusion_baselines = {35,15};

    // One can also use the map for exclusion baseline:
    // r.map_exclusion_threshold_baselines = {{"Ch0.", 10}, {"Ch7.", 11}, {"Ch10.", 11}, {"Ch17.", 9}, {"Ch20.", 9}, {"Ch27.", 12}, {"Ch30.", 10}, {"Ch37.", 10},};
    r.exclusion_window = 500; // time in ns that it will jump for baseline
    r.filter = 9; // denoise filter.. if filter = 0, no denoise is done.
    r.OnlySomeEvents = false; // Do you want only a few events? Choose it wisely (you can set stopEvent)
    r.stopEvent = 1000; //
    r.noBaseline = false; // if you dont want to calculate baseline
    // r.saveFilter = true;

    // r.channels = {1};
    // r.channels = {1,2,3,4};
      

    
    r.adc_read_all_data();

    
}

