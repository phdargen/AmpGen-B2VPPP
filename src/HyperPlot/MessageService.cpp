#include "AmpGen/HyperPlot/MessageService.h"


MessageSerivce* MessageSerivce::s_messageService = 0;

///Static function to retrive the MessageSerivce singleton.
///At the same time, also change the state of the MessageSerivce.
MessageSerivce& MessageSerivce::getMessageService(ErrorType errorType){
  
  //if the singleton doesn't exist, return make it
  if (s_messageService == 0){
  	s_messageService = new MessageSerivce();
  }
  
  //See if we should be printing the error type given.
  //Set the _printOrNot variable accordingly 
  s_messageService->_printOrNot = s_messageService->_outputOptions[errorType];
  


  //Set the error type to the one given
  if ( s_messageService->_printOrNot){

    //if the message type has changed, and no endl was called, then 
    //force a new line.
    if (s_messageService->_errorType != errorType && s_messageService->_endlCalled == 0){
      *s_messageService << std::endl;
    }

  	s_messageService->_errorType = errorType;
  }
  
  //return the singleton
  return *s_messageService;

}
 
///Static function to retrive the MessageSerivce singleton
///without altering its message state 
MessageSerivce& MessageSerivce::getMessageService(){
  
  //if the singleton doesn't exist, return make it
  if (s_messageService == 0){
  	s_messageService = new MessageSerivce();
  }
  
  //return the singleton
  return *s_messageService;

}


///constuctor
///
MessageSerivce::MessageSerivce() :
  _stream(std::cout),
  _endlCalled(true),
  _errorCount(0)
{

  _outputOptions[WELCOME] = false;
  _outputOptions[ERROR  ] = true ;
  _outputOptions[INFO   ] = true ;
  _outputOptions[VERBOSE] = false;
  _outputOptions[GOODBYE] = false;

  _outputHeaders[WELCOME] = "HyperPlot Welcome : ";
  _outputHeaders[ERROR  ] = "\033[31mHyperPlot Error : ";
  _outputHeaders[INFO   ] = "\033[32mHyperPlot Info : ";
  _outputHeaders[VERBOSE] = "HyperPlot Verbose : ";
  _outputHeaders[GOODBYE] = "HyperPlot Goodbye : ";

}

///print how many errors messages have been sent to the output stream
///
void MessageSerivce::printErrorCount(){

  //turns out that it never gets destructed so this doesn't work. Bit of a shame
  if (_errorCount != 0){
  	_errorCount--; //This isn't actually an error, so -1
  	getMessageService(ERROR) << "There were " << _errorCount << " errors during runtime" << std::endl;
  }	
}

///destuctor
///
MessageSerivce::~MessageSerivce(){
  
  //turns out that it never gets destructed so this doesn't work. Bit of a shame
  if (_errorCount != 0){
    _errorCount--;
  	getMessageService(ERROR) << "There was " << _errorCount << " errors during runtime" << std::endl;
  }

}


