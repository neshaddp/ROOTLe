//////////////////////////////////////////////////////////
//
//    Evt2Cal.cc
//    Created       : 2014/05/22 16:53:56
//--------------------------------------------------------
//    Comment : Convert raw data (NSCL DAQ 10.2-104)
//              to calibrated/calculated ROOT TTree
//--------------------------------------------------------
//
//--------------------------------------------------------
//    Modified to parse data in format nscldaq 11.0
//    J. Pereira, March 2016
// ------------------------------------------------------------------------------------------------------
//    The new (current) ddas data format has THREE importnat differences
//    with respect to previous versions:
//
//    1. Besides the fragment header + ring-item header + body included in 
//    previous nscldaq versions, there is a BODY HEADER 
//   
//    2. Unlike in previous versions, the time-stamp included in the fragment and 
//    body headers is calibrated (in ns units).
//
//    3. The body includes a new 32-bit word with a "module-identifying" word
//    defined in the Readout function CMyEventSegment::Read (CMyEventSegment.cpp)
// -------------------------------------------------------------------------------------------------------
//   
//   RING ITEM from EVENT BUILDER (ringbuffer cedaqb) from LENDA DDAS:                                     
//   (see http://docs.nscl.msu.edu/daq/newsite/nscldaq-11.0/x4493.html) 
//   ("NEW" corresponds to parts of the ring-item not present in versions prior to 11.0) 
//   -----------------------------------------------------------------------------------------------------
//   | RING-ITEM | 1 Self-inclusive size payload         (1 32-bit word) 
//   | HEADER    | 2 Ring-item type                      (1 32-bit word)
//   |----------------------------------------------------------------------------------------------------
//   | BODY      | 1 Size of header (=20 bytes)          (1 32 bit word) 
//   | HEADER    | 2 Time stamp                          (2 32-bit words)
//   | (NEW)     | 3 Source ID                           (1 32 bit word)
//   |           | 4 Barrier type                        (1 32 bit word)
//   |----------------------------------------------------------------------------------------------------
//   | BODY      | 1 Self-inclusive total body size in bytes  (1 32 bit word)
//   |           | ---------------------------------------------------------------------------------------
//   |           | FRAGMENT 0    | FRAGMENT  | 1 Time stamp           (2 32-bit words)
//   |           | (from e.g     | HEADER    | 2 Source ID            (1 32 bit word)
//   |           | channel 0     |           | 3 Payload size         (1 32 bit word)
//   |           | of LENDA      |           | 4 Barrier type         (1 32 bit word)  
//   |           | DDAS)         |------------------------------------------------------------------------
//   |           |               | PAYLOAD   | RING-ITEM | 1 Self-inclusive size payload (1 32-bit word)
//   |           |               | FRAGMENT  | HEADER    | 2 Ring-item type              (1 32-bit word) 
//   |           |               |           |------------------------------------------------------------
//   |           |               |           | BODY      | 1 Size of header (=20 bytes)  (1 32 bit word) 
//   |           |               |           | HEADER    | 2 Time stamp                  (2 32-bit words)
//   |           |               |           | (NEW)     | 3 Source ID                   (1 32 bit word)
//   |           |               |           |           | 4 Barrier type                (1 32 bit word)  
//   |           |               |           |------------------------------------------------------------
//   |           |               |           | BODY      | (This is where event data are located)
//   |           |               |           |           | "NEW" Module identifying word (1 32 bit word) 
//   |           |               |           |           | Rest of data from DDAS channel....
//   |           |----------------------------------------------------------------------------------------
//   |           | FRAGMENT 1    | FRAGMENT  | 1 Time stamp           (2 32-bit words)
//   |           | (from e.g     | HEADER    | 2 Source ID            (1 32 bit word)
//   |           | channel 1     |           | 3 Payload size         (1 32 bit word)
//   |           | of LENDA      |           | 4 Barrier type         (1 32 bit word)  
//   |           | DDAS)         |------------------------------------------------------------------------
//   |           |               | PAYLOAD   | RING-ITEM | 1 Self-inclusive size payload  (1 32-bit word)
//   |           |               | FRAGMENT  | HEADER    | 2 Ring-item type               (1 32-bit word) 
//   |           |               |           |-------------------------------------------------------------
//   |           |               |           | BODY      | 1 Size of header (=20 bytes)   (1 32 bit word) 
//   |           |               |           | HEADER    | 2 Time stamp                   (2 32-bit words)
//   |           |               |           | (NEW)     | 3 Source ID                    (1 32 bit word)
//   |           |               |           |           | 4 Barrier type                 (1 32 bit word)  
//   |           |               |           |-------------------------------------------------------------
//   |           |               |           | BODY      | (This is where event data are located)
//   |           |               |           |           | "NEW" Module identifying word (1 32 bit word) 
//   |           |               |           |           | Rest of data from DDAS channel....
//   |           |----------------------------------------------------------------------------------------
//   |           | .nth-FRAGMENT |  ........................................................
//   |----------------------------------------------------------------------------------------------------
// 
//   The program checks the first word after the ring-item header. If the value of that word is 20, then
//   it assumes that there is a body header (the first word in a body header is the size of the header, which should
//   be equal to 20 bytes), and parses the rest of the ring-item accordingly. If the value is different than 20, then
//   the program assumes that there is no body header, and treats the read variable as the total payload size 
//   (as it would do for data formats generated by nscldaq versions previous to 11.0). 
//   The idea is that the new program knows how to parse data from versions 11.0 and older.





#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <sys/stat.h>
#include <cstdlib>

#include "TFile.h"
#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "S800Event.hh"
#include "DDASEvent.hh"
#include "LendaPacker.hh"
#include "R00TLeSettings.hh"

#include "S800Settings.hh"
#include "S800Calibration.hh"

#include "DAQdefs.h"

#include "Utilities.hh"

using namespace std;

int main(int argc, char* argv[])
{
   if (argc != 3  && argc !=4) {
      Error("Evt2Cal","Usage: Evt2Cal InputFile OutputFile [RunNumber]");
      return 0;
   }
   
   double time_start = get_time();  
   TStopwatch* sw;
   sw = new TStopwatch();
   sw->Start(kTRUE);
   signal(SIGINT,signalhandler);

   Int_t RunNumber=-1;
   if (argc == 4 ){ // A RunNumber was given
     RunNumber = atoi(argv[3]);
   }

   // Get the name of the input file from the arguments
   TString InputFile  = TString(argv[1]);

   FileStat_t info;
   if (gSystem->GetPathInfo(InputFile, info) != 0) {
      Error("Evt2Cal", "Input file %s does not exist.", InputFile.Data());
      return 0;
   }

   // Open the input file.
   FILE *infile;
   infile = fopen(InputFile.Data(),"r");
   
   Info("Evt2Cal", "Input file:  %s", InputFile.Data());
   Info("Evt2Cal", "Input file size: %7.2f MB", (float)info.fSize/1024./1024.);

   // Get the name of the output file from the arguments
   TString OutputFile = TString(argv[2]);
   
   // Open the output file.
   TFile *outfile = new TFile(OutputFile,"RECREATE");
   if (outfile->IsZombie()) {
      return 0;
   }

   Info("Evt2Cal", "Output file: %s", OutputFile.Data());

   // Initialize the data structures for the event building.
   // Needs cleanup
   Int_t   buffers = 0;
   Int_t   buffers_physev = 0;
   int64_t bytes_read = 0;
   size_t  count_read;

   int8_t   buffer08[4096];
   uint16_t buffer16[4096];
   uint32_t buffer32[4096];

   // for S800
   uint16_t evtlength;
   uint16_t *pevent;
   S800* s800 = new S800;

   // for ddaschannel
   uint32_t *body_ptr;
   ddaschannel* dchan = new ddaschannel;

   // for NSCLDAQ-format data
   size_t count;
   uint32_t size_of_entire_item;
   uint32_t type_of_entire_item;

   uint32_t size_of_built_body_header;
   uint64_t timestamp_built_body_header;
   uint32_t sourceid_built_body_header;
   uint32_t barrier_built_body_header;

   uint32_t total_size_of_body;
   uint32_t bytes_read_in_body;
   uint64_t timestamp;
   uint32_t sourceid;
   uint32_t fragment_payload_size;
   uint32_t barrier;
   uint32_t size_of_ring_item;
   uint32_t type_of_ring_item;


   uint32_t size_of_body_header;
   uint64_t timestamp_body_header;
   uint32_t sourceid_body_header;
   uint32_t barrier_body_header;

   uint32_t dummy32;
   uint16_t dummy16;

   bool display = false;;
   bool bodyhead = false;

   uint64_t s800_tstamp;
   // uint64_t ddas_tstamp;

   // prepare a tree instead of DecodedEvent
   int nentries;
   //TTree* outtree   = new TTree("rawtree","S800 and DDAS raw events");
   TTree* outtree   = new TTree("caltree","S800 and DDAS calibrated events");
   S800Event* s800event = new S800Event;
   DDASEvent* ddasevent = new DDASEvent;
   //outtree->Branch("s800event", &s800event, 320000);
   //outtree->Branch("ddasevent", &ddasevent, 320000);


   // Here is where the object TheSettings is created
   R00TLeSettings * TheR00TLeSettings = new R00TLeSettings();

   // S800Calc branch
   S800Calc*  s800calc    = new S800Calc;
   outtree->Branch("s800calc",   &s800calc,   320000);

   // LendaEvent branch
   LendaPacker *thePacker = new LendaPacker(TheR00TLeSettings);
   // thePacker->SetFilter(2,0,2,1);//FL FG d w
   // thePacker->SetObjectFilter(2,0,2,4);//Fl fg d w
   thePacker->SetGates(25,11,25,11);
   thePacker->SetTraceDelay(120);
   thePacker->FindAndSetMapAndCorrectionsFileNames(RunNumber);

   LendaEvent* lendaevent = new LendaEvent;
   outtree->Branch("lendaevent", &lendaevent, 320000);

   outtree->BranchRef();

   // Parameters and settings
   TString prmdirname = gSystem->Getenv("R00TLe_PRM");
   TString prmfilename("Evt2Cal.prm");
   TString prmfilenamefull = prmdirname + "/" + prmfilename;
   std::ifstream prmfile(prmfilenamefull);
   Info("Evt2Cal", "Reading a parameter file %s", prmfilenamefull.Data());
   if (prmfile.fail()) {
      Error("Evt2Cal", "No such prm file %s", prmfilenamefull.Data());
      return 4;
   }

   S800Settings    *set = new S800Settings(prmfilenamefull.Data());
   S800Calibration *cal = new S800Calibration(set);

   // initialize
   nentries    = 0;
   s800event->Clear();
   ddasevent->Reset();
   s800calc  ->Clear();
   lendaevent->Clear();

   Bool_t BadLastEvent=kFALSE;
   //Loop over the entirety of the input file.
   while (!signal_received) {
     //first try and read something from the file
     //if it is at the end then you must try and read 
     //from the file inorder for the eof bit to true
     
     // The header of built main ring item (from the event builder): 2 uint32_t's.
     // A 32-bit size of entire item (in bytes)
     count_read  = fread(&size_of_entire_item, sizeof(uint32_t), 1, infile);   
     if ( feof(infile) ){ // if the end of tile bit has been set
       if ( bytes_read != info.fSize){
	 Error("Evt2Cal","Total number of bytes read %ld not equal to file size %lld at the end of file parsing",bytes_read,info.fSize);
       }
       break; //end main unpacking loop
     }
     bytes_read += 1 * sizeof(uint32_t);

     if (display) cout << " SIZE (dec) of EVB ring item (RI) " << dec << size_of_entire_item << dec << 
     ". Bytes read: " << bytes_read << endl;
     
        


     // A 32-bit type of item 
     count_read  = fread(&type_of_entire_item, sizeof(uint32_t), 1, infile);
     bytes_read += 1 * sizeof(uint32_t);
     if (display) cout << " TYPE of EVB RI  " << hex << type_of_entire_item << dec << 
     ". Bytes read: " << bytes_read << endl;
     
     
     // The pointer should be pointing at the beginning of the main "body" built ring item
     // Switching according to the type of the entire event.
     // Doing nothing except for physics events.
     switch (type_of_entire_item) {
        case PHYSICS_EVENT:

	  if (display) cout << "-------- We got a physics event -----" << endl;
	  
	  bytes_read_in_body = 0; // bytes_read_in_body is used to check that all data in ring-item body is read 

	  buffers_physev++; //Count new physics-event buffer
	  

	  // In nscldaq-11, the next uint32_t (dummy32) is the size of the main ring-item body header in bytes (it should be 20)
	  // In nscldaq-10, the next uing32_t (dummy32) is the total body size from the body 
	  count_read  = fread(&dummy32, sizeof(uint32_t), 1, infile);
	  bytes_read +=  sizeof(uint32_t);
	  if (display) cout << "SIZE of EVB RI body-header (=20) or TOTAL SIZE of BODY (if nscldaq10): " << dec << dummy32 << dec << 
			 ". Bytes read: " << bytes_read << endl;

	  // If we are in the body header, the dummy word should give the header-size, which,
	  // for a PhYSICS_EVENT type must be 20 bytes.
	  // We want to verify it to make sure that there are indeed body headers. If the size is not 20,
	  // we'll assume that we are dealing with data format prior to nscldaq-11.0
	  // J. Pereira, March 2016
	  if (dummy32 == 20) {
	    bodyhead = true;
	  } else {
	    bodyhead = false;
	  }

	  if (bodyhead) { // If there is body header, read it before the total payload size

	    if (display) cout << "Parsing body header of main ring item: " << endl;
	  
	    size_of_built_body_header = dummy32; // dummy word is the size of the body header

	    // 64-bit timestamp
	    count_read  = fread(&timestamp_built_body_header,sizeof(uint64_t), 1, infile);
	    bytes_read += sizeof(uint64_t);
	    if (display) cout << "TIMESTAMP from body header of EVB RI   " << hex << timestamp_built_body_header << dec << 
			   ". Bytes read: " << bytes_read << endl;
	    
	    // 32-bit sourceid
	    count_read  = fread(&sourceid_built_body_header,sizeof(uint32_t), 1, infile);
	    bytes_read += sizeof(uint32_t);
	    if (display) cout << "SOURCE ID from body header of EVB RI   " << hex << sourceid_built_body_header << dec << 
			   ". Bytes read: " << bytes_read << endl;
	    
	    // 32-bit barrier 
	    count_read  = fread(&barrier_built_body_header,sizeof(uint32_t), 1, infile);
	    bytes_read += sizeof(uint32_t);
	    if (display) cout << "BARRIER from body header of EVB RI   " << hex << barrier_built_body_header << dec << 
			   ". Bytes read: " << bytes_read << endl;
	       
	    // Total (selfinclusive) payload size (in bytes) of entire ring-item body from event builder 
	    count_read  = fread(&total_size_of_body, sizeof(uint32_t), 1, infile);
	    bytes_read += sizeof(uint32_t);
	    bytes_read_in_body += sizeof(uint32_t);
	    if (display) cout << "TOTAL body SIZEfrom body of EVB RI (dec) " << dec << total_size_of_body << dec << 
			   ". Bytes read: " << bytes_read << endl;
 
	  } else { // If there is no body-header, then the dummy word is the total payload size of the entire ring-item body 
	    total_size_of_body = dummy32; 
	    bytes_read_in_body += sizeof(uint32_t);
	    bytes_read = bytes_read - 1*sizeof(uint32_t); //total body size is self-inclusive, but we counted it when we read dummy32
	  }

	  if (display) cout << "Bytes read: " << dec << bytes_read <<
	    "  Total body size: " << total_size_of_body <<
	    "  File size: " << info.fSize << endl;
	  

	  if ((bytes_read+total_size_of_body) > info.fSize){
	    BadLastEvent=kTRUE;
	    cout<<"\n\nThis fragment goes over the total files size"<<endl;
	    cout<<"Skipping and ending the main loop"<<endl;
	    break;//Break out of Switch over Type of event
	  }
	  
	  
	  while(kTRUE) {

  	       // Reading a fragment ///////////////////////////////////
	    if (display) cout << "---------------------We are now reading fragments --- " << endl; 

	       // 64-bit timestamp
	       count_read  = fread(&timestamp,sizeof(uint64_t), 1, infile);
	       bytes_read += sizeof(uint64_t);

	       if (display) cout << "TIMESTAMP from FRAGMENT HEADER   " << hex << timestamp << dec << 
			      ". Bytes read: " << bytes_read << endl;
	      

	       // 32-bit sourceid
	       count_read  = fread(&sourceid,sizeof(uint32_t), 1, infile);
	       bytes_read += sizeof(uint32_t);

	       if (display) cout << "SOURCE ID from FRAGMENT HEADER  " << hex << sourceid << dec << 
			      ". Bytes read: " << bytes_read << endl;
	      

	       // 32-bit fragment payload size (in bytes)
	       count_read  = fread(&fragment_payload_size,sizeof(uint32_t), 1, infile);
	       bytes_read += sizeof(uint32_t);

	       if (display) cout << "PAYLOAD SIZE from FRAGMENT HEADER (dec) " << dec << fragment_payload_size  << dec << 
			      ". Bytes read: " << bytes_read << endl;
	      


	       // 32-bit barrier 
	       count_read  = fread(&barrier,sizeof(uint32_t), 1, infile);
	       bytes_read += sizeof(uint32_t);

	       if (display) cout << "BARRIER from body FRAGMENT HEADER  " << hex << barrier << dec << 
			      ". Bytes read: " << bytes_read << endl;

	       

	       // Reading ring-item header in fragment ////////////////////////

	       // 32-bit total size of ring item (in bytes)
	       count_read  = fread(&size_of_ring_item,sizeof(uint32_t), 1, infile);
	       bytes_read += sizeof(uint32_t);
	       if (display) cout << "SIZE from FRAGMENT RI HEADER  " << dec << size_of_ring_item << dec << 
			      ". Bytes read: " << bytes_read << endl;
   
       

	       // 32-bit type of ring item
	       count_read  = fread(&type_of_ring_item,sizeof(uint32_t), 1, infile);
	       bytes_read += sizeof(uint32_t);
	       if (display) cout << "TYPE from FRAGMENT RI HEADER  " << hex << type_of_ring_item << dec << 
			      ". Bytes read: " << bytes_read << endl;




	       // Decode the body of event differently depending on which type of data is identified in the header.
	       // For each, pass the data from the file into the DecodedEvent to be read.

	       switch (sourceid) {
		  case DDAS_SOURCEID:
		     switch (type_of_ring_item) {
			case PHYSICS_EVENT:
			{

			  if (display) cout << "-------- Now parsing fragment body of DDAS event -----" << endl;
			  
   			  if (bodyhead) { // If there is a body header, read it
			  
			    
			    // Size of the body header in bytes (it should be 20)
			    count_read  = fread(&size_of_body_header, sizeof(uint32_t), 1, infile);
			    bytes_read +=  sizeof(uint32_t);
			    if (display) cout << "SIZE from FRAGMENT BODY HEADER  " << dec << size_of_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;
			    		    
			    // 64-bit timestamp
			    count_read  = fread(&timestamp_body_header,sizeof(uint64_t), 1, infile);
			    bytes_read += sizeof(uint64_t);
			    if (display) cout << "TIMESTAMP from FRAGMENT BODY HEADER  " << hex << timestamp_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;
			    
			    // 32-bit sourceid
			    count_read  = fread(&sourceid_body_header,sizeof(uint32_t), 1, infile);
			    bytes_read += sizeof(uint32_t);
			    if (display) cout << "SOURCE ID from FRAGMENT BODY HEADER  " << hex << sourceid_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;
			    
			    // 32-bit barrier 
			    count_read  = fread(&barrier_body_header,sizeof(uint32_t), 1, infile);
			    bytes_read += sizeof(uint32_t);
			    if (display) cout << "BARRIER from FRAGMENT BODY HEADER  " << hex << barrier_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;

			    
			    // 7 sizeof(uint32_t)'s need to be skipped from the beginning of fragment ring item to get to the body
			    // [7 = 2 (ring-item header) +  5 (body header: 3 32-bit words + 1 64-bit word)]
			    count       = size_of_ring_item - 7*sizeof(uint32_t);
			    
			    
			  } else {
			    // If there is no body header, just skip 2 sizeof(uint32_t)'s from the ring-item header
			    count       = size_of_ring_item - 2*sizeof(uint32_t);
			    
			  }
			  

			   count_read  = fread(buffer32, sizeof(int8_t), count, infile);
			   bytes_read += count_read * sizeof(int8_t);


			   // The pointer should now be pointing at the end of the ring item.
			   // buffer32[0] should have the address of the beginning of the DDAS event body
			   body_ptr = buffer32;
			 
			   //dchan->Reset();
			   dchan = new ddaschannel; // Needs to be careful with push_back with class object...
			   dchan->UnpackChannelData(body_ptr);

                           // std::cout << dchan->GetTimeLow() + dchan->GetTimeHigh() * TMath::Power(2,32) << "\t";
                           // std::cout << (int)(s800_tstamp - ddas_tstamp)
                           //           << std::endl;

			   // Add decoded ddaschannel instance to fDDASEvent
			   ddasevent->AddChannelData(dchan);
			   // (ddasevent->GetData()).emplace_back(dchan); // This may also work...

			   break;

			}
			default:
			{
			   // In case of non-physics events, just translate the pointer
			   // to the end of this ring item. 

			   if (bodyhead) { // If there is bodyheader

			     // 7 sizeof(uint32_t)'s need to be skipped from the beginning of ring item to get to the body
			     // [7 = 2 (ring-item header) +  5 (body header: 3 32-bit words + 1 64-bit word)]
			     count       = size_of_ring_item - 7*sizeof(uint32_t);
			     bytes_read += 5 * sizeof(uint32_t); // Count 5 32-bit words in body header. JP
			     
			   } else {
			     // If there is no body header, just skip 2 sizeof(uint32_t)'s from the ring-item header
			     count       = size_of_ring_item - 2*sizeof(uint32_t);
			     
			   }

			  count_read  = fread(buffer32, sizeof(int8_t), count, infile);
			  bytes_read += count_read * sizeof(int8_t);
	
	
			}
		     }
		     break;
	  
		  case S800_SOURCEID:
		  {

		    switch (type_of_ring_item) {
			case PHYSICS_EVENT:
			{
			  
			  if (display) cout << "-------- Now parsing fragment body of S800 event -----" << endl;

			  
   			  if (bodyhead) { // If there is a body header, read it
			  
			    // Size of the body header in bytes (it should be 20)
			    count_read  = fread(&size_of_body_header, sizeof(uint32_t), 1, infile);
			    bytes_read +=  sizeof(uint32_t);
			    if (display) cout << "SIZE from FRAGMENT BODY HEADER  " << dec << size_of_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;
			    		    
			    // 64-bit timestamp
			    count_read  = fread(&timestamp_body_header,sizeof(uint64_t), 1, infile);
			    bytes_read += sizeof(uint64_t);
			    if (display) cout << "TIMESTAMP from FRAGMENT BODY HEADER  " << hex << timestamp_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;
			    
			    // 32-bit sourceid
			    count_read  = fread(&sourceid_body_header,sizeof(uint32_t), 1, infile);
			    bytes_read += sizeof(uint32_t);
			    if (display) cout << "SOURCE ID from FRAGMENT BODY HEADER  " << hex << sourceid_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;	    

			    // 32-bit barrier 
			    count_read  = fread(&barrier_body_header,sizeof(uint32_t), 1, infile);
			    bytes_read += sizeof(uint32_t);
			    if (display) cout << "BARRIER from FRAGMENT BODY HEADER  " << hex << barrier_body_header << dec << 
					   ". Bytes read: " << bytes_read << endl;
			    
			  }
  
			  		  
			  // Store s800_tstamp in case this fragment is from S800.
			  s800_tstamp = timestamp_body_header;
			  if (display) cout << "Time stamp: " << hex << s800_tstamp << endl;
			   
			  // Check now selfinclusive body payload (16-bit event length)
			  count_read = fread(&evtlength, sizeof(uint16_t), 1, infile);
			  if (count_read == 0) {
			    Warning("Evt2Cal","S800 header but no data.");
			    break;
			  }
			  bytes_read += sizeof(uint16_t);


			  
			  evtlength--; // selfinclusive

			  if (display) cout << "evtlength:  " << evtlength << endl;
			   
			   // (evtlength-1) 16-bit words for the data

			  
			   count_read = fread(buffer16, sizeof(uint16_t), evtlength, infile);
			   if (count_read == 0) {
			      Warning("Evt2Cal","S800 expected %u words but found none.",evtlength);
			      break;
			   }
			   bytes_read += evtlength*sizeof(uint16_t);

			   // The pointer should now be pointing at the end of the ring item.
			   // buffer16[0] should have the address of the beginning of the S800 data to be analyzed


			   
			   pevent = buffer16;

			   s800->Clear();
			   int error = s800->DecodeS800(pevent,evtlength);

			   if(error){
			      Warning("Evt2Cal", 
				      "An error (%d) occurred in DecodeS800() while processing %s. Continuing...",
				      error,
				      InputFile.Data());
			      continue;
			   }
			   
			   // Add decoded S800 instance to fS800
			   s800event->SetS800(*s800);
			   break;
			}
			default:
			{
			   // In case of non-physics events, just translate the pointer
			   // to the end of this ring item



			  if (bodyhead) { // If there is bodyheader
			    
			    // 7 sizeof(uint32_t)'s need to be skipped from the beginning of ring item to get to the body
			    // [7 = 2 (ring-item header) +  5 (body header: 3 32-bit words + 1 64-bit word)]
			    count       = size_of_ring_item - 7*sizeof(uint32_t);
			    bytes_read += 5 * sizeof(uint32_t); // Count 5 32-bit words in body header. JP
			    
			  } else {
			    // If there is no body header, just skip 2 sizeof(uint32_t)'s from the ring-item header
			    count       = size_of_ring_item - 2*sizeof(uint32_t);
			    
			  }
			  
			  count_read  = fread(buffer08, sizeof(int8_t), count, infile);
			  bytes_read += count_read * sizeof(int8_t);
			  
			  
			}
		     }//end switch type of ring item
		  }//end case S800SourceId
		  break;
	  
		  // If source is neither DDAS nor S800, just skip the data
		  default:
		     count_read = fread(buffer08, sizeof(int8_t),
				   size_of_ring_item - 2 * sizeof(uint32_t),
				   infile);
		     bytes_read += (size_of_ring_item-2*sizeof(uint32_t))*sizeof(int8_t);
	       }

	       // Data read in bytes from the beginning of fragments
	       bytes_read_in_body +=     sizeof(uint64_t) // size of the holder for timestamp
		                     +   sizeof(uint32_t) //                        sourceid
		                     +   sizeof(uint32_t) //                        fragment payload size
		                     +   sizeof(uint32_t) //                        barrier
		                     +   fragment_payload_size; // actual fragment payload size
	
	       // check if all the event fragments are processed
	       if (total_size_of_body == bytes_read_in_body) {
		  // Build the calibrated S800 object from data in the uncalibrated S800 object
		  // cal->S800Calculate(s800event->GetS800(),s800calc); // this should work
		  s800calc->ApplyCalibration(s800event->GetS800(),cal); // this should work as well
                  // Setting the timestamp

                  s800calc->SetTS((long long int)s800_tstamp);
		  //long long tsreturn = s800calc->GetTS();
		  //cout << "Set TS S800: " << s800_tstamp << "  Get TS S800: " << tsreturn << endl; 
		  
		    
                  // Just for check...
                  /*
		  for (uint i = 0; i < ddasevent->GetData().size(); i++) {
                    dchan = ddasevent->GetData()[i];

                    std::cout << "tstamp_diff:\t"
                              << (int)(s800calc->GetTS()
                                       - (dchan->GetTimeLow() + dchan->GetTimeHigh() * TMath::Power(2,32)))
                              << std::endl;
                  }
                  */
		  
		  thePacker->MakeLendaEvent(lendaevent, ddasevent,nentries);
		  
		  lendaevent->Finalize();
		  
		  // populate tree; clear event for the next;
		  outtree->Fill();
		  // reinitialize
		  s800event->Clear();

                  // Avoiding memory leak
		  for (uint i = 0; i < ddasevent->GetData().size(); i++) {
		    delete ddasevent->GetData()[i];
		  }
		  ddasevent->GetData().clear(); //Clearing Vector
		  
		  s800calc  ->Clear();
		  lendaevent->Clear();
		  nentries++;

		  break;
	       } else {
		  // Reading the next fragment...
	       }
	    } // Done reading all the fragments.
	    break;
    
	 default:
	    // anything other than PHYSICS_EVENT (state change, scalers, etc.)
	    count       = size_of_entire_item - 2*sizeof(uint32_t);
	    count_read  = fread(buffer08, sizeof(int8_t), count, infile);
	    // The pointer should now be pointing at the beginning of the next built ring item.
	    // Increment the total bytes read.
	    bytes_read += count_read * sizeof(int8_t);
      }

      // Write the trees out to disk every denom events.
      const Int_t denom = 10000;
      buffers++;
      if(buffers % denom == 0){
	outtree->AutoSave();
     
	double time_end = get_time();
	Progress("Evt2Cal", "%8d buffers %5.0f MB (%6.2f%%) read at %5.2f buffers/sec... %6.2f sec to go\r",
		 buffers, 
		 float(bytes_read/(1024*1024)), 
		 float(bytes_read/(1024*1024)) / float(info.fSize/1024./1024.) * 100.,
		 float(buffers/(time_end-time_start)),
		 ((float)(info.fSize-bytes_read)/(float)bytes_read)*(time_end-time_start)
		 );
      }
      
      if (BadLastEvent==kTRUE){
	break;
      }

   }

   Info("Evt2Cal","Total of %d data buffers (%5.2f MB read)",
	buffers,
	(float)(bytes_read/(1024.*1024.)));

   Info("Evt2Cal","Total of %d physics_event buffers written",
	buffers_physev);

   Info("Evt2Cal","Total of %d raw events written",
	nentries);

   double time_end = get_time();
   Info("Evt2Cal","Decoded %f buffers/sec.",
	(float)(buffers/(time_end - time_start)));

   sw->Stop();
   sw->Print("u");
  
   //Final cleanup and writing of files.
   outtree->Write("",TObject::kOverwrite);
   TheR00TLeSettings->Write();
   outfile->Close();
   // If outfile cannot be closed because the file size exceeds fgMaxTreeSize, try TTree::SetMaxTreeSize

   return 0;
}
