/*****************************************/
/* This channel initializes a new        */
/* DDAS channel event.                   */
/*****************************************/

#include "ddaschannel.hh"


ClassImp(ddaschannel);

ddaschannel::ddaschannel() : TObject() {
  channelnum = 0;
  chanid = 0;
  slotid = 0;
  crateid = 0;
  channelheaderlength = 0;
  channellength = 0;
  finishcode = 0;
  overflowcode = 0;
  tracelength = 0;
  timelow = 0;
  timehigh = 0;
  timecfd = 0;
  cfdtrigsourcebit=0;  //added CP
  energy = 0;
  
  time = 0;
  cfd = 0;
  
  has_esum = false;
  has_qdc = false;
  has_exclock = false;
  

  exclocklow = 0; // Include external clock. JP 1/21/2016
  exclockhigh = 0;
  exclock = 0;

  energySums.reserve(4);
  qdcSums.reserve(8);
  trace.reserve(200);
}

ddaschannel::ddaschannel(const ddaschannel& obj)
  : TObject(obj), 
    channelnum(obj.channelnum),
    chanid(obj.chanid),
    slotid(obj.slotid),
    crateid(obj.crateid),
    channelheaderlength(obj.channelheaderlength),
    channellength(obj.channellength),
    finishcode(obj.finishcode),
    overflowcode(obj.overflowcode),
    tracelength(obj.tracelength),
    timelow(obj.timelow),
    timehigh(obj.timehigh),
    timecfd(obj.timecfd),
    cfdtrigsourcebit(obj.cfdtrigsourcebit),  //added CP
    energy(obj.energy),
    time(obj.time),
    cfd(obj.cfd),
    has_esum(obj.has_esum),
    has_qdc(obj.has_qdc),
    has_exclock(obj.has_exclock),
    exclocklow(obj.exclocklow), // Include external clock. JP 1/21/2016 
    exclockhigh(obj.exclockhigh),
    exclock(obj.exclock),
    energySums(obj.energySums),
    qdcSums(obj.qdcSums),
    trace(obj.trace)
{} 

ddaschannel& ddaschannel::operator=(const ddaschannel& obj)
{
  if (this!=&obj) {
    channelnum = obj.channelnum;
    chanid = obj.chanid;
    slotid = obj.slotid;
    crateid = obj.crateid;
    channelheaderlength = obj.channelheaderlength;
    channellength = obj.channellength;
    finishcode = obj.finishcode;
    overflowcode = obj.overflowcode;
    tracelength = obj.tracelength;
    timelow = obj.timelow;
    timehigh = obj.timehigh;
    timecfd = obj.timecfd;
    cfdtrigsourcebit = obj.cfdtrigsourcebit;  //added CP
    energy = obj.energy;
    time = obj.time;
    cfd = obj.cfd;
    has_esum = obj.has_esum;
    has_qdc = obj.has_qdc;
    has_exclock = obj.has_exclock;
    exclocklow = obj.exclocklow; // Include external clock. JP 1/21/2016 
    exclockhigh = obj.exclockhigh;
    exclock = obj.exclock;
    energySums = obj.energySums;
    qdcSums = obj.qdcSums;
    trace = obj.trace;
  }

  return *this;
}

void ddaschannel::Reset() {
  channelnum = 0;
  chanid = 0;
  slotid = 0;
  crateid = 0;
  channelheaderlength = 0;
  channellength = 0;
  finishcode = 0;
  tracelength = 0;
  cfdtrigsourcebit = 0;  //Added CP
  
  timelow = 0;
  timehigh = 0;
  timecfd = 0;
  energy = 0;
  
  time = 0;
  cfd = 0;

  has_esum = false;
  has_qdc = false;
  has_exclock = false;
 
  exclocklow = 0; // Include external clock. JP 1/21/2016
  exclockhigh = 0;
  exclock = 0;

  energySums.clear();
  qdcSums.clear();
  trace.clear();
}

ddaschannel::~ddaschannel() {


}

void ddaschannel::SetChannelID(UInt_t data){
  chanid = ((Int_t) ((data & CHANNELIDMASK)));
}

void ddaschannel::SetSlotID(UInt_t data){
  slotid = ((Int_t) ((data & SLOTIDMASK) >> 4));
}

void ddaschannel::SetCrateID(UInt_t data){
  crateid = ((Int_t) ((data & CRATEIDMASK) >> 8));
}

void ddaschannel::SetChannelHeaderLength(UInt_t data){
  channelheaderlength = ((Int_t) ((data & HEADERLENGTHMASK) >> 12 ));
}

void ddaschannel::SetChannelLength(UInt_t data){
  channellength = ((Int_t) ((data & CHANNELLENGTHMASK) >> 17));
}

void ddaschannel::SetOverflowCode(UInt_t data){
  overflowcode = ((Int_t) ((data & OVERFLOWMASK) >> 30));
}

void ddaschannel::SetFinishCode(UInt_t data){
  finishcode = ((Int_t) ((data & FINISHCODEMASK) >> 31 ));
}

void ddaschannel::SetID(UInt_t data){
  //broken function
}

void ddaschannel::SetTimeLow(UInt_t data){
  timelow = data;
}

void ddaschannel::SetTimeHigh(UInt_t data){
  timehigh = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddaschannel::SetTimeCFD(UInt_t data){
  timecfd = ((Int_t) ((data & BIT29to16MASK) >> 16));
}

void ddaschannel::SetCFDTriggerSourceBit(UInt_t data){
  cfdtrigsourcebit = ((Int_t) ((data & BIT30MASK) >> 30 ));  //modified CP
}

void ddaschannel::SetTime(){
  
  time = timecfd/16384.0 + 2*(timelow + timehigh * 4294967296.0) - cfdtrigsourcebit;
}

void ddaschannel::SetEnergy(UInt_t data){
  energy = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddaschannel::SetTraceLength(UInt_t data){
  //tracelength = ((Int_t) ((data & UPPER16BITMASK) >> 16));
  
  // If there is an overflow, bit 31 of the word tracelength + energy is set to 1. 
  tracelength = ((Int_t) ((data & UPPER15BITMASK) >> 16));
  Int_t overflow = ((Int_t) ((data & UPPER16BITMASK) >> 31));
  //if (overflow == 1) cout << " ***** Warning: detector overflow " << endl; 
		      
 }

void ddaschannel::SetEnergySums(UInt_t data){
  energySums.push_back(data);
}

void ddaschannel::SetQDCSums(UInt_t data){
  qdcSums.push_back(data);
}



void ddaschannel::SetExClockLow(UInt_t data){
  exclocklow = data;
}

void ddaschannel::SetExClockHigh(UInt_t data){
  exclockhigh = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddaschannel::SetExClock(){
  exclock = exclocklow + exclockhigh * 4294967296.0;
  //cout << exclock << endl;
}



void ddaschannel::SetTraceValues(UInt_t data){
  //cout << "set trace values " << data << endl;
  trace.push_back((data & LOWER16BITMASK));
  trace.push_back((data & UPPER16BITMASK)>>16);
  //cout << "set trace values1 " << data << endl;
}

void ddaschannel::UnpackChannelData(const uint32_t *data){

  //put all variables in known state before unpacking
  Reset();

  //first four words of data identifiers, times, and energies
  //cout << "data ----------" << endl;
  //cout << dec << "data0   : " << data[0] << endl;
  //cout << hex << "data1   : " << data[1] << endl;
  //cout << hex << "data2   : " << data[2] << endl;
  //cout << hex << "data3   : " << data[3] << endl;
  //cout << hex << "data4   : " << data[4] << endl;
  //cout << hex << "data5   : " << data[5] << endl;

  
  //data[0] is the number of half words in the event, skip it
  *data++;


  // A module-identifying 32-bit word is included in current ddas version (March 2016)
  // We'll skip it for now. J. Pereira
  uint32_t module_identifier = *data;
  //cout << "Module identifier: " << hex << module_identifier << endl;
  //if (module_identifier == MODIDENTIFIER) {
  //cout << "Module identifier: " << hex << module_identifier << endl;
  data++; // First word is Module-identifying word -- skip it. JP, March 2016
  //} 



  //Using the first word of DDAS information extract channel identifiers
  SetChannelID(*data);
  SetSlotID(*data);
  SetCrateID(*data);
  SetChannelHeaderLength(*data);
  SetChannelLength(*data);
  SetOverflowCode(*data);
  SetFinishCode(*data++);

  //Second word
  // set most significant bits of time stamp
  SetTimeLow(*data++);

  //Third word
  // set least significant bits of time stamp
  SetTimeHigh(*data);

  //set cfd time
  SetTimeCFD(*data);

  //set CFDTriggerSourceBit
  SetCFDTriggerSourceBit(*data++); //added CP

  //set the full time
  SetTime();

  //Fourth word
  //set energy
  SetEnergy(*data);

  //set trace length
  SetTraceLength(*data++);

  // finished upacking the minimum set of data


  // Determine if there are additional header words to unpack
  if(channelheaderlength == 4) {
    has_exclock = false;  
    has_esum = false;    
    has_qdc = false;     

  } else if(channelheaderlength == 6) {
    has_exclock = true;  // External clock
    has_esum = false;    
    has_qdc = false;  

  } else if (channelheaderlength == 8) { 
    has_exclock = false;  
    has_esum = true;    // Energy sum and baselines
    has_qdc = false;     

  } else if (channelheaderlength == 12) {
    has_exclock = false;  
    has_esum = false;    
    has_qdc = true;     // QDC sums

  } else if (channelheaderlength == 10) {
    has_exclock = true;  // External clock
    has_esum = true;    // Energy sum and baselines
    has_qdc = false;     

  } else if (channelheaderlength == 14) {
    has_exclock = true;  // External clock
    has_esum = false;    
    has_qdc = true;     // QDC sums

  } else if (channelheaderlength == 16) {
    has_exclock = false;  
    has_esum = true;    // Energy sum and baselines
    has_qdc = true;     // QDC sums

  } else {
    cout << "Inconsistent head lengths  " << channelheaderlength  << endl;
    cout << "This is a fatal error - data is corrupted, please investigate " << endl;
    exit(1);
  }

  if (has_esum) { // Unpack four words with energy sums and baselines
    for(int z=0; z<4; z++){
      SetEnergySums(*data++);
    }
  }

  if (has_qdc) { // Unpack eight words with QDC sums
    for(int z=0; z<8; z++){
      SetQDCSums(*data++);
    }
  }

  if (has_exclock) { // Unpack two words with external clock
    
    //cout << "External clock detected" << endl;
    
    SetExClockLow(*data++);
    SetExClockHigh(*data++);
    SetExClock();
  }







  //determine if more data should be unpacked based on a comparison between channel header length
  //and channel length
  if(channelheaderlength != channellength){
    //cout << "extra info unpack " <<endl;
    //more unpacking data
    if(channellength != (channelheaderlength + tracelength/2)){
      cout << "Inconsistent lengths between channel length - " << channellength
           << " , header length - " << channelheaderlength 
           << " , and trace length - " << tracelength << endl;
      cout << "This is a fatal error - data is corrupted, please investigate " << endl;
      exit(1);
    }

    // if trace length is non zero, retrieve the trace
    if(tracelength !=0) {
      //cout << "trace unpack "<<endl;
      for(Int_t z = 0; z < tracelength/2; z++){
	//upack two 12 bit trace values from 32 bit data word
	SetTraceValues(*data++);
      }   
    } //finished unpacking trace
    
  } //finished unpacking extra data


}
