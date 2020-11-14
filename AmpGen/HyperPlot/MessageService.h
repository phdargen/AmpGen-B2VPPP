/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * This is used to print all of the output from HyperPlot
 * to a chosen output stream (by default std::cout).
 * Access to this class is done through the precompiler
 * definitions.
 *
 *  WELCOME_LOG  - A welcome message in the constructor of each class
 *  ERROR_LOG    - Use to print error messages to the screen 
 *  INFO_LOG     - Use to print useful information to the screen
 *  VERBOSE_LOG  - Use for verbose information
 *  GOODBYE_LOG  - A goodbye message in the destructor of each class
 *
 *  ERROR_COUNT  - Print out how many error messages have been shown.
 *                 This is useful to call at the end of the main function
 * 
 *
 **/
 
#ifndef MESSAGE_SERVICE_HH
#define MESSAGE_SERVICE_HH

#include <iostream>
#include <sstream>
#include <map>

#include "TString.h"


class MessageSerivce{
  
  private:
  
  ///static variable that holds the singleton object
  static MessageSerivce* s_messageService;
  
  MessageSerivce();
  
  public:
  
  ///Define the possible error types
  enum ErrorType{WELCOME, ERROR, INFO, VERBOSE, GOODBYE};

  ///the output stream - by default this is std:cout
  std::ostream& _stream;

  ///this map decides what message types should be printed 
  ///to the output stream (verbose, welcome and goodbye
  ///are off by default)
  std::map<ErrorType, bool>    _outputOptions;

  ///this is the string that proceeds each of the 
  ///message types
  std::map<ErrorType, TString> _outputHeaders;
  
  ///bool was the last command sent to the steam
  ///a std::cout ?
  bool _endlCalled;

  ///For the current message type, should I be printing
  ///to the stream (means you don't have to use the map
  ///many times in a row i.e. faster)
  bool _printOrNot;

  ///The current message type of the MessageSerivce 
  ///
  ErrorType _errorType;
  
  ///Count the number of errors that happened.
  ///
  long int _errorCount;

  static MessageSerivce& getMessageService(ErrorType errorType);
  static MessageSerivce& getMessageService();

  void printErrorCount();
  
  ///This makes is possible the do messageSerivce << 
  ///
  template <class T>MessageSerivce  &operator<< (const T &v) 
  { 
  
    if (_printOrNot) {
      if (_endlCalled){
        _stream << _outputHeaders[_errorType];
        if (_errorType == ERROR) _errorCount++;
        _endlCalled = 0;
      }
      _stream << v;
    }
    return *this;
  }
  
  //This allows you to make a custom endl (MessageSerivce::endl).
  //Don't think I need it, but interesting

  // function that takes a custom stream, and returns it
  //typedef MessageSerivce& (*MessageSerivceManipulator)(MessageSerivce&);
  

  // take in a function with the custom signature
  // MessageSerivce& operator<<(MessageSerivceManipulator manip)
  // {
  //   call the function, and return it's value
  //   return manip(*this);
  //  }

  // define the custom endl for this stream.
  // note how it matches the `MessageSerivceManipulator`
  // function signature
  //  static MessageSerivce& endl(MessageSerivce& stream)
  //  {
  // print a new line
  //   std::cout << std::endl;

  // do other stuff with the stream
  // std::cout, for example, will flush the stream
  //   stream << "Called MessageSerivce::endl!" << std::endl;

  //   return stream;
  //  }

  //The following makes endl work


  /// this is the type of std::cout
  ///
  typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

  /// this is the function signature of std::endl
  ///
  typedef CoutType& (*StandardEndLine)(CoutType&);

  /// define an operator<< to take in std::endl
  ///
  MessageSerivce& operator<<(StandardEndLine manip)
  { 
    //SAM -> I'm assuming that when this gets called it's 
    //always endl - I imagine this isn't always the case.
    if (_printOrNot){
      _endlCalled = true;

    // call the function, but we cannot return it's value
      _stream << "\033[0m";
      manip(_stream);
      
    }
    return *this;
  }

  ~MessageSerivce(); 

};

#define WELCOME_LOG MessageSerivce::getMessageService(MessageSerivce::WELCOME)
#define ERROR_LOG   MessageSerivce::getMessageService(MessageSerivce::ERROR  )
#define INFO_LOG    MessageSerivce::getMessageService(MessageSerivce::INFO   )
#define VERBOSE_LOG MessageSerivce::getMessageService(MessageSerivce::VERBOSE)
#define GOODBYE_LOG MessageSerivce::getMessageService(MessageSerivce::GOODBYE)

#define ERROR_COUNT MessageSerivce::getMessageService().printErrorCount();

#endif
