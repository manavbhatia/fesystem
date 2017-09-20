// $Id: TimeLogs.h,v 1.1.4.4 2007-05-08 05:19:14 manav Exp $

#ifndef __fesystem_time_logs_h__
#define __fesystem_time_logs_h__

// C/C++ includes
#include <map>
#include <stack>
#include <string>
#include <memory>
#include <ostream>
#include <sys/time.h>

#ifdef HAVE_LOCALE
#include <locale>
#endif


// FESystem includes


// libMesh includes
#include "utils/perf_log.h"
#include "utils/o_string_stream.h"



namespace FESystemUtility
{
  /** 
  * This class has been adapted from the libMesh \p FESystemUtility::PerformanceLogging class to change 
  * some behavior.
  * The \p FESystemUtility::PerformanceLogging class allows monitoring of specific events.
  * An event is defined by a unique string that functions as
  * a label.  Each time the event is executed data are recorded.
  * This class is particulary useful for finding performance
  * bottlenecks. 
  *
  */
  
  // ------------------------------------------------------------
  // FESystemUtility::PerformanceLogging class definition
  class PerformanceLogging
  {
    
public:
    
    /**
    * Constructor.  \p label_name is the name of the object, which
     * will bw print in the log to distinguish it from other objects.
     * \p log_events is a flag to optionally
     * disable logging.  You can use this flag to turn off
     * logging without touching any other code.
     */
    PerformanceLogging(const std::string& label_name="",
                       const bool log_events=true);
    
    /**
    * Destructor. Calls \p clear() and \p print_log().
     */
    ~PerformanceLogging();
    
    /**
      * Clears all the internal data and returns the
     * data structures to a pristine state.  This function
     * checks to see if it is currently monitoring any
     * events, and if so errors.  Be sure you are not
     * logging any events when you call this function.
     */
    void clear();
    
    /**
      * Start monitoring the event named \p label.
     */
    void start_event(const std::string &label,
                     const std::string &header="");
    
    /**
      * Stop monitoring the event named \p label.
     */
    void stop_event(const std::string &label,
                    const std::string &header="");
    
    /**
      * Suspend monitoring of the event. 
     */
    void pause_event(const std::string &label,
                     const std::string &header="");
    
    /**
      * Restart monitoring the event.
     */
    void restart_event(const std::string &label,
                       const std::string &header="");
    
    /**
      * @returns a string containing:
     * (1) Basic machine information (if first call)
     * (2) The performance log
     */
    std::string get_log() const;
    
    /**
      * @returns a string containing ONLY the information header.
     */
    std::string get_info_header() const;
    
    /**
      * @returns a string containing ONLY the log information
     */
    std::string get_perf_info() const;
    
    /**
      * Print the log.
     */
    void print_log() const;
    
//    /**
//      * @returns the total time spent on this event.
//     */
//    double get_total_time() const
//      {return total_time;}
    
    
private:
      
      
      /**
      * The label for this object.
       */
      const std::string label_name;
    
    /**
      * Flag to optionally disable all logging.
     */
    const bool log_events;
    
    /**
      * The time we were constructed or last cleared.
     */
    struct timeval tstart;
    
    /**
      * The actual log.
     */
    std::map<std::pair<std::string,
      std::string>,
      PerfData> log;
    
    
    /**
      * The log of time spent in each header
     */
    std::map<std::string, PerfData> header_time;

    
    /**
      * Flag indicating if print_log() has been called.
     * This is used to print a header with machine-specific
     * data the first time that print_log() is called.
     */
    static bool called;
    
    /**
      * Prints a line of 'n' repeated characters 'c'
     * to the output string stream "out".
     */
    void _character_line(const unsigned int n,
                         const char c,
                         OStringStream& out) const;
  };
  
  
  
  // ------------------------------------------------------------
  // FESystemUtility::PerformanceLogging class inline member funcions
  
  // Typedefs we might need
#ifdef HAVE_LOCALE
  typedef std::ostreambuf_iterator<char, std::char_traits<char> > TimeIter;
  typedef std::time_put<char, TimeIter> TimePut;
#endif
  
  
  
  
  /// this class stores the information for a header
  class HeaderEventData
    {
public: 
      /// enumeration of the status for event logs
      enum EventStatus {RUNNING, STOPPED};
      
      /// constructor
      HeaderEventData();
      
      /// destructor
      ~HeaderEventData();
      
      /// @returns current event
      const std::string& getCurrentEvent() const;
      
      /// @returns call count of an event
      unsigned int getEventCallCount(const std::string& event_name) const;
      
      /// @returns current status
      FESystemUtility::HeaderEventData::EventStatus getCurrentStatus() const;
      
      /// removes the current event
      void removeCurrentEvent(const std::string& event_name);
      
      /// sets the current event
      void setCurrentEvent(const std::string& event_name);
      
      /// @returns true if some events exist in this header
      bool eventsExist() const;
      
      /// prints the event stack for this header
      void printEventCallCount(std::ostream& output);
      
protected:
        /// current event
        std::stack<std::string> event_stack;
      
      /// call count for open events
      std::map<std::string, unsigned int> event_call_count;
    };
  
  /// this class builds on top of the libmesh FESystemUtility::PerformanceLogging class, and provides some 
  /// additional functionality
  class TimeLogs
    {
public:
      /// constructor
      TimeLogs();
      
      /// destructor
      ~TimeLogs();
      
      /// this method sets the event under the given header. This behaves in the following way:
      /// 
      /// -- if no other event is being logged in the given header, then this event 
      /// defines the first event and is started
      ///
      /// -- if the current event in the header being logged is the same as the given event,
      /// then it checks if the event was paused or not. If it was paused, it restarts the 
      /// event. Otherwise, it adds to the call counter of the event
      /// 
      /// -- if the current event in the header is different, then if the event is not already 
      /// paused, it pauses the event, and sets the current logged event to the given event
      /// 
      void setEvent(const std::string& event_name, const std::string& header);
      
      
      /// this method unsets the current event and switches the event back to the one that was 
      /// being logged before the current event. It is an error to unset an event if it is not 
      /// currently being logged.
      void unsetEvent(const std::string& event_name, const std::string& header);
      
      
protected:
        
        /// the performance log object
        std::auto_ptr<FESystemUtility::PerformanceLogging> performance_logging;
      
      /// map of headers and events
      std::map<std::string, FESystemUtility::HeaderEventData*> event_map;
    };
  
}



#endif // __fesystem_time_logs_h__
